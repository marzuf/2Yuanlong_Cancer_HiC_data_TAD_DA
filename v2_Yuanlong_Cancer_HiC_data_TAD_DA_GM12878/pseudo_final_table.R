#!/usr/bin/Rscript

outFolder <- file.path("PSEUDO_FINAL_TABLE")
dir.create(outFolder, recursive=TRUE)

options(scipen=100)

startTime <- Sys.time()

# Rscript pseudo_final_table.R

plotType <- "svg"

source("settings.R")


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480 , 7)
myWidth <- ifelse(plotType == "png", 600, 10)
plotCex <- 1.4

require(flux)
require(foreach)
require(doMC)
registerDoMC(nCpu)

pipFolder <- "."
pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")

cat(paste0("... pipFolder = ", pipFolder, "\n"))

minTADsize <- 3

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script11same_name <- "11sameNbr_runEmpPvalCombined"

setDir <- ""
entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
entrez2symb_dt$symbol <- as.character(entrez2symb_dt$symbol)
stopifnot(!duplicated(entrez2symb_dt$entrezID))

all_result_dt <- foreach(hicds = all_obs_hicds, .combine='rbind') %dopar% {
  exprds = all_obs_exprds[[paste0(hicds)]][1]
  hicds_dt <- foreach(exprds = all_obs_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    cat(paste0("... start ", hicds, " -  ", exprds, "\n"))
    
    ### RETRIEVE THE GENE2TAD ASSIGNMENT
    g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
    
    cat(paste0("g2tFile = ", g2tFile, "\n"))
    
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    
    
    ### RETRIEVE THE TAD ASSIGNMENT
    tadFile <- file.path(pipFolder, hicds, "genes2tad", "all_assigned_regions.txt")
    stopifnot(file.exists(tadFile))
    tad_DT <- read.delim(tadFile, header=F, col.names = c("chromo","region", "start", "end"), stringsAsFactors = FALSE)
    tad_DT <- tad_DT[grepl("_TAD", tad_DT$region),]
    stopifnot(nrow(tad_DT) > 0)
    
    ### RETRIEVE THE GENES USED IN THE PIPELINE - script0
    stopifnot(dir.exists(file.path(pipOutFolder, hicds, exprds, script0_name)))
    geneListFile <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneListFile))
    pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted
    stopifnot(pipeline_geneList %in% g2t_DT$entrezID)    
    g2t_DT <- g2t_DT[g2t_DT$entrezID %in% pipeline_geneList,]
    
    fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
    stopifnot(file.exists(fc_file))
    tad_fc <- eval(parse(text = load(fc_file)))
    all_regs <- names(tad_fc)
    
    stopifnot(all_regs %in% g2t_DT$region)
    stopifnot(all_regs %in% tad_DT$region)
    all_regs_start <- setNames(tad_DT$start, tad_DT$region)
    all_regs_end <- setNames(tad_DT$end, tad_DT$region)
    
    corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
    stopifnot(file.exists(corr_file))  
    tad_corr <- eval(parse(text = load(corr_file)))
    stopifnot(setequal(all_regs, names(tad_corr)))
    
    comb_empPval_file <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
    stopifnot(file.exists(comb_empPval_file))
    comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
    stopifnot(setequal(all_regs, names(comb_empPval)))
    comb_empPval <- comb_empPval[all_regs]
    # ADJUST THE PVAL
    tad_adjCombPval <- p.adjust(comb_empPval, method="BH")
    stopifnot(names(tad_adjCombPval) == all_regs)
    
    tad_genes <- sapply(all_regs, function(x) {
      x_genes <- g2t_DT$entrezID[g2t_DT$region == x]
      stopifnot(length(x_genes) >= minTADsize)
      stopifnot(x_genes %in% entrez2symb_dt$entrezID)
      x_symbols <- entrez2symb_dt$symbol[entrez2symb_dt$entrezID %in% x_genes]
      stopifnot(length(x_symbols) == length(x_genes))
      x_symbols <- sort(x_symbols)
      paste0(x_symbols, collapse=",")
    })
    names(tad_genes) <- all_regs
    
    stopifnot(all_regs %in% names(all_regs_end))
    stopifnot(all_regs %in% names(all_regs_start))
    
    
    exprds_hicds_dt <- data.frame(
      hicds = hicds,
      exprds = exprds,
      region = all_regs,
      
      start = all_regs_start[all_regs],
      end = all_regs_end[all_regs],
      
      region_genes = tad_genes[all_regs],
      meanLogFC = tad_fc[all_regs],
      meanCorr = tad_corr[all_regs],
      
      adjPvalComb = tad_adjCombPval[all_regs],
      
      stringsAsFactors = FALSE
    )
    rownames(exprds_hicds_dt) <- NULL
    
    
    
    exprds_hicds_dt <- exprds_hicds_dt[order(exprds_hicds_dt$adjPvalComb),]
    
    exprds_hicds_dt
  } # end-foreach iterating over exprds
  hicds_dt
} # end-foreach iterating over hicds
outFile <- file.path(outFolder, "all_result_dt.Rdata")
save(all_result_dt, file = outFile)
cat(paste0("... written: ", outFile,  "\n"))

















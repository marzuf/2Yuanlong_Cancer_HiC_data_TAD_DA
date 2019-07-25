### CREATE TABLE
# -----------------------
#  TAD_ID | genes | meanFC | meanCorr | adj. comb. p-value | signif. at FDR 0.1 ?  | signif. at FDR 0.2  

# Rscript create_final_table.R

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "create_final_table.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script19_name <- "19_SAM_emp_measurement"
script8c_name <- "8c10000_runAllDown"

require(foreach)
require(doMC)
registerDoMC(40)

pipFolder<- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "CREATE_FINAL_TABLE"
dir.create(outFolder, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
hicds <- args[1]
exprds <- args[2]

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
entrez2symb_dt$symbol <- as.character(entrez2symb_dt$symbol)
stopifnot(!duplicated(entrez2symb_dt$entrezID))


fdr_thresh1 <- 0.1
fdr_thresh2 <- 0.2
corr_type <- "sample_meanCorr_allDS"

buildTable <- TRUE

samp_type <- "sameNbr"
fixKbSize <- ""

corrEmpFDR_folder <- file.path("SAM_EMP_MEASUREMENT_MEANCORR", samp_type, fixKbSize)
stopifnot(dir.exists(corrEmpFDR_folder))

notAdjCombPval_folder <- file.path("CREATE_EMPPVALCOMB", samp_type, fixKbSize)
stopifnot(dir.exists(notAdjCombPval_folder))


if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}

minTADsize <- 3


hicds = all_hicds[1]

if(buildTable) {
  
  
  all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %do% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start building DT for :", hicds, " - ", exprds, "\n")
      
      
      ### RETRIEVE THE GENE2TAD ASSIGNMENT
      g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      
      ### RETRIEVE THE GENES USED IN THE PIPELINE - script0
      script0_name <- "0_prepGeneData"
      stopifnot(dir.exists(file.path(pipOutFolder, hicds, exprds, script0_name)))
      geneListFile <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneListFile))
      pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted
      stopifnot(pipeline_geneList %in% g2t_DT$entrezID)    
      g2t_DT <- g2t_DT[g2t_DT$entrezID %in% pipeline_geneList,]
      
      corrEmpFDR_file <- file.path(corrEmpFDR_folder, hicds, exprds, "all_empFDR.Rdata")
      stopifnot(file.exists(corrEmpFDR_file))
      corrEmpFDR <- eval(parse(text = load(corrEmpFDR_file)))
      
      fcEmpFDR_file <- file.path(pipOutFolder, hicds, exprds, script19_name, "empFDR_list.Rdata")
      stopifnot(file.exists(fcEmpFDR_file))
      fc_empFDR <- eval(parse(text = load(fcEmpFDR_file)))
      
      fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
      stopifnot(file.exists(fc_file))
      tad_fc <- eval(parse(text = load(fc_file)))
      all_regs <- names(tad_fc)
      
      corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
      stopifnot(file.exists(corr_file))  
      tad_corr <- eval(parse(text = load(corr_file)))
      
      stopifnot(setequal(all_regs, names(tad_corr)))
      
      notAdjCombPval_file <- file.path(notAdjCombPval_folder, hicds, exprds, "meanCorr_meanLogFC_notAdjCombEmpPval.Rdata")
      stopifnot(file.exists(notAdjCombPval_file))
      tad_notAdjCombPval <- eval(parse(text = load(notAdjCombPval_file)))
      stopifnot(setequal(all_regs, names(tad_notAdjCombPval)))
      tad_adjCombPval <- p.adjust(tad_notAdjCombPval, method="BH")
      
      fc_FDR <- fc_empFDR[["empFDR_logFC"]]
      # the 1st FC value for which FDR is smaller than the thresh
      fc_FDR_co1 <- min(as.numeric(as.character(na.omit(names(fc_FDR)[fc_FDR <= fdr_thresh1]))))
      # will return Inf is none below the thresh
      fc_FDR_co2 <- min(as.numeric(as.character(na.omit(names(fc_FDR)[fc_FDR <= fdr_thresh2]))))
      
      corr_FDR <- corrEmpFDR[[paste0(corr_type)]][["empFDR"]]
      corr_FDR_co1 <- min(as.numeric(as.character(na.omit(names(corr_FDR)[corr_FDR <= fdr_thresh1])))) # will return Inf is none below the thresh
      corr_FDR_co2 <- min(as.numeric(as.character(na.omit(names(corr_FDR)[corr_FDR <= fdr_thresh2])))) # will return Inf is none below the thresh
      
      FDR_signif_co1 <- setNames(abs(tad_fc[all_regs]) >= fc_FDR_co1 &
                                   tad_corr[all_regs] >= corr_FDR_co1, all_regs)
      
      stopifnot(setequal(all_regs, names(FDR_signif_co1)))
      
      
      FDR_signif_co2 <- setNames(abs(tad_fc[all_regs]) >= fc_FDR_co2 &
                                   tad_corr[all_regs] >= corr_FDR_co2, all_regs)
      
      
      
      stopifnot(setequal(all_regs, names(FDR_signif_co2)))
      
      
      # RETRIEVE RATIO DOWN
      rd_file <- file.path(pipOutFolder, hicds, exprds, script8c_name, "all_obs_ratioDown.Rdata")
      stopifnot(file.exists(rd_file))
      tad_rD <- eval(parse(text = load(rd_file)))
      stopifnot(setequal(names(tad_rD), all_regs))
      
      
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
      
      
      exprds_hicds_dt <- data.frame(
        hicds = hicds,
        exprds = exprds,
        region = all_regs,
        region_genes = tad_genes[all_regs],
        meanLogFC = tad_fc[all_regs],
        meanCorr = tad_corr[all_regs],
        
        ratioDown = tad_rD[all_regs],
        
        meanLogFC_FDR_co1 = fc_FDR_co1,
        meanLogFC_FDR_co2 = fc_FDR_co2,
        
        meanCorr_FDR_co1 = corr_FDR_co1,
        meanCorr_FDR_co2 = corr_FDR_co2,
        
        
        adjPvalComb = tad_adjCombPval[all_regs],
        signifFDRco1 = FDR_signif_co1[all_regs],
        signifFDRco2 = FDR_signif_co2[all_regs],
        stringsAsFactors = FALSE
      )
      rownames(exprds_hicds_dt) <- NULL
      colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "signifFDRco1"] <- paste0("signifFDR_",
                                                                                       fdr_thresh1)
      colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "signifFDRco2"] <- paste0("signifFDR_",
                                                                                       fdr_thresh2)
      
      colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "meanLogFC_FDR_co1"] <- paste0("meanLogFC_thresh_FDR",
                                                                                            fdr_thresh1)
      colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "meanLogFC_FDR_co2"] <- paste0("meanLogFC_thresh_FDR",
                                                                                            fdr_thresh2)
      
      
      colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "meanCorr_FDR_co1"] <- paste0("meanCorr_thresh_FDR",
                                                                                            fdr_thresh1)
      colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "meanCorr_FDR_co2"] <- paste0("meanCorr_thresh_FDR",
                                                                                            fdr_thresh2)
      
      
      exprds_hicds_dt <- exprds_hicds_dt[order(exprds_hicds_dt$adjPvalComb),]
      
      exprds_hicds_dt
    } # end-foreach iterating over exprds
    hicds_dt
  } # end-foreach iterating over hicds
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  save(all_result_dt, file = outFile)
  cat(paste0("... written: ", outFile,  "\n"))
  
  
  takecols <- colnames(all_result_dt)[! colnames(all_result_dt) %in% c(paste0("meanLogFC_thresh_FDR",
                                                                              fdr_thresh1),
                                                                       paste0("meanLogFC_thresh_FDR",
                                                                              fdr_thresh2),
                                                                       paste0("meanCorr_thresh_FDR",
                                                                              fdr_thresh1),
                                                                       paste0("meanCorr_thresh_FDR",
                                                                              fdr_thresh2))]
  
  
  all_result_dt_txt <- all_result_dt[, takecols]
  outFile <- file.path(outFolder, "all_result_dt.txt")
  write.table(all_result_dt_txt, file = outFile, sep="\t", quote=F, row.names=F, col.names=T,append=F)
  cat(paste0("... written: ", outFile,  "\n"))

    
} else {
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  all_result_dt <- eval(parse(text = load(outFile)))
  
  

  
  
  
}



#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))










options(scipen=100)

# Rscript boundary_effect_gene.R   # to run all datasets in one shot

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "boundary_effect_gene.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

tiesMeth <- "min"

buildData <- TRUE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4


pipFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

corrPvalsFolder <- file.path("SAMPLE_MEANCORR_EMPPVALS")
stopifnot(dir.exists(corrPvalsFolder))

outFolder <- "BOUNDARY_EFFECT_GENE"
dir.create(outFolder, recursive=TRUE)

corr_type <- "meanCorr"
samp_type <- "sameNbr"
permut_type <- "allDS"
fixSizeKb <- ""
sampName <- paste0("empPval-", samp_type, "-", corr_type, " - ", permut_type)



pvalFolder <- "CREATE_EMPPVALCOMB"
samp_type <- "sameNbr"
fixSizeKb <- ""

stopifnot(dir.exists(pvalFolder))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2 | length(args) == 0)
hicds <- args[1]
exprds <- args[2]

if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}


entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
entrez2symb_dt$strand <- as.character(entrez2symb_dt$strand)
stopifnot(entrez2symb_dt$strand == "+" | entrez2symb_dt$strand == "-")



if(buildData) {
  
  ### BUILD SIGNIF ALONG FDR THRESH
  cat("... start building signif. along FDR thresh\n")
  all_gene_tad_expr_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds_gene_tad_expr_DT <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      # RETRIEVE PVAL COMBINED DATA
      
      pvalFile <- file.path(pvalFolder, samp_type, fixSizeKb, hicds, exprds,
                            "meanCorr_meanLogFC_notAdjCombEmpPval.Rdata")
      stopifnot(file.exists(pvalFile))
      combPvals <- eval(parse(text = load(pvalFile)))
      adj_combPvals <- p.adjust(combPvals, method="BH")
      adj_combPvals_rank <- rank(adj_combPvals, ties = tiesMeth) # best rank = smallest value
      stopifnot(adj_combPvals[names(adj_combPvals_rank)[adj_combPvals_rank == 1]] == min(adj_combPvals))
      stopifnot(names(adj_combPvals) == names(adj_combPvals_rank))
      
      tad_rank_DT <- data.frame(region=names(adj_combPvals_rank),
                                TAD_rank = adj_combPvals_rank, 
                                TAD_adjCombPval = adj_combPvals,                          
                                stringsAsFactors = FALSE)
      
      
      de_file <- file.path(pipOutFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")
      stopifnot(file.exists(de_file))
      de_DT <- eval(parse(text = load(de_file)))
      de_DT$genes <- as.character(de_DT$genes)
      g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      geneList_file <- file.path(pipOutFolder, hicds, exprds,script0_name,   "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      geneList <- eval(parse(text = load(geneList_file)))
      stopifnot(names(geneList) %in% de_DT$genes)
      stopifnot(geneList %in% g2t_DT$entrezID)
      de_DT <- de_DT[de_DT$genes %in% names(geneList),]
      stopifnot(nrow(de_DT) == length(geneList))
      de_DT$entrezID <- sapply(de_DT$genes, function(x) as.character(geneList[as.character(x)]))
      stopifnot(!is.na(de_DT$entrezID))
      
      gene_DE <- setNames(de_DT$adj.P.Val, de_DT$entrezID)
      gene_logFC <- setNames(de_DT$logFC, de_DT$entrezID)
      gene_DE_rank <- rank(gene_DE, ties = tiesMeth) # best rank = smallest value
      stopifnot(gene_DE[names(gene_DE_rank)[gene_DE_rank == 1]] == min(gene_DE))
      stopifnot(names(gene_DE_rank) == names(gene_logFC))
      stopifnot(names(gene_DE_rank) == names(gene_DE))
      gene_rank_DT <- data.frame(entrezID=names(gene_DE_rank), 
                                 gene_rank = gene_DE_rank, 
                                 gene_adjPval = gene_DE,
                                 gene_logFC = gene_logFC,
                                 stringsAsFactors = FALSE)
      
      g2tFile <- file.path(pipFolder, hicds,  "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "gene_start", "gene_end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      stopifnot(geneList %in% g2t_DT$entrezID)
      g2t_DT <- g2t_DT[g2t_DT$entrezID %in% geneList,]
      
      stopifnot(setequal(g2t_DT$entrezID, gene_rank_DT$entrezID))
      g2t_rank_DT <- merge(gene_rank_DT, g2t_DT[,c("entrezID", "region", "gene_start", "gene_end")], by="entrezID")
      
      stopifnot(setequal(tad_rank_DT$region, g2t_rank_DT$region))
      
      gene_tad_expr_DT <- merge(tad_rank_DT, g2t_rank_DT, by ="region")
      
      gene_tad_expr_DT$hicds <- hicds
      gene_tad_expr_DT$exprds <- exprds
      
      
      
      ### RETRIEVE THE TAD POSITIONS
      tadposFile <- file.path(pipFolder,hicds, "genes2tad", "all_assigned_regions.txt")
      stopifnot(file.exists(tadposFile))
      tadpos_DT <- read.delim(tadposFile, header=F, col.names=c("chromo", "region", "region_start", "region_end"), stringsAsFactors = FALSE)
      stopifnot(is.numeric(tadpos_DT$region_start))
      stopifnot(is.numeric(tadpos_DT$region_end))
      tadpos_DT <- tadpos_DT[grepl("_TAD", tadpos_DT$region),,drop=FALSE] 
      tadpos_DT <- tadpos_DT[tadpos_DT$region %in% g2t_DT$region,]
      tadpos_DT$chromo <- NULL
      
      gene_tad_expr_DT_temp <- merge(gene_tad_expr_DT, tadpos_DT, by ="region")
      
      
      gene_tad_expr_DT <- merge(gene_tad_expr_DT_temp, entrez2symb_dt[,c("entrezID", "strand")], by = "entrezID", all.x=TRUE)
      stopifnot(!is.na(gene_tad_expr_DT))
      stopifnot(nrow(gene_tad_expr_DT_temp) == nrow(gene_tad_expr_DT))
      
      
      gene_tad_expr_DT$tss_start <- ifelse(gene_tad_expr_DT$strand == "+", gene_tad_expr_DT$gene_start, 
                                           ifelse(gene_tad_expr_DT$strand == "-", gene_tad_expr_DT$gene_end, NA))
      stopifnot(!is.na(gene_tad_expr_DT$tss_start))
      
      gene_tad_expr_DT$tss_end <- ifelse(gene_tad_expr_DT$strand == "+", gene_tad_expr_DT$gene_end, 
                                           ifelse(gene_tad_expr_DT$strand == "-", gene_tad_expr_DT$gene_start, NA))
      stopifnot(!is.na(gene_tad_expr_DT$tss_end))
      
      
      stopifnot(is.numeric(gene_tad_expr_DT$gene_start))
      stopifnot(is.numeric(gene_tad_expr_DT$gene_end))
      stopifnot(is.numeric(gene_tad_expr_DT$region_start))
      stopifnot(is.numeric(gene_tad_expr_DT$region_end))
      
      
      
      
      
      
      gene_tad_expr_DT$gene_midPos <- (gene_tad_expr_DT$gene_start+gene_tad_expr_DT$gene_end)/2
      
      gene_tad_expr_DT$region_midPos <- (gene_tad_expr_DT$region_start+gene_tad_expr_DT$region_end)/2
      
      gene_tad_expr_DT$minDistMidToBD <- pmin( abs(gene_tad_expr_DT$gene_midPos-gene_tad_expr_DT$region_start),
                                            abs(gene_tad_expr_DT$gene_midPos-gene_tad_expr_DT$region_end))
      
      
      gene_tad_expr_DT$distMidToMid <- abs(gene_tad_expr_DT$gene_midPos-gene_tad_expr_DT$region_midPos)



      gene_tad_expr_DT$minDistStartToBD <- pmin( abs(gene_tad_expr_DT$tss_start-gene_tad_expr_DT$region_start),
                                            abs(gene_tad_expr_DT$tss_start-gene_tad_expr_DT$region_end))
      
      
      gene_tad_expr_DT$distStartToMid <- abs(gene_tad_expr_DT$tss_start-gene_tad_expr_DT$region_midPos)


      gene_tad_expr_DT$minDistEndToBD <- pmin( abs(gene_tad_expr_DT$tss_end-gene_tad_expr_DT$region_start),
                                            abs(gene_tad_expr_DT$tss_end-gene_tad_expr_DT$region_end))
      
      
      gene_tad_expr_DT$distEndToMid <- abs(gene_tad_expr_DT$tss_end-gene_tad_expr_DT$region_midPos)
      
      
      gene_tad_expr_DT
      
    } # end-foreach iterating over exprds
    exprds_gene_tad_expr_DT
  } # end-foreach iterating over hicds
  
  outFile <- file.path(outFolder, "all_gene_tad_expr_dt.Rdata")  
  save(all_gene_tad_expr_dt, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else { # end-if build data
  outFile <- file.path(outFolder, "all_gene_tad_expr_dt.Rdata")  
  all_gene_tad_expr_dt <- eval(parse(text = load(outFile)))
}

all_gene_tad_expr_dt$subtype_col <- all_cols[paste0(all_cmps[paste0(all_gene_tad_expr_dt$exprds)])]

all_gene_tad_expr_dt$dataset <- paste0( all_gene_tad_expr_dt$hicds, " - ", all_gene_tad_expr_dt$exprds)
all_datasets <- unique(all_gene_tad_expr_dt$dataset) 
nDS <- length(unique( all_datasets ))

all_gene_tad_expr_dt$minDistMidToBD_log10 <- log10(all_gene_tad_expr_dt$minDistMidToBD)
all_gene_tad_expr_dt$distMidToMid_log10 <- log10(all_gene_tad_expr_dt$distMidToMid)


all_gene_tad_expr_dt$minDistStartToBD_log10 <- log10(all_gene_tad_expr_dt$minDistStartToBD)
all_gene_tad_expr_dt$distStartToMid_log10 <- log10(all_gene_tad_expr_dt$distStartToMid)

all_gene_tad_expr_dt$minDistEndToBD_log10 <- log10(all_gene_tad_expr_dt$minDistEndToBD)
all_gene_tad_expr_dt$distEndToMid_log10 <- log10(all_gene_tad_expr_dt$distEndToMid)


###################################
################################### PLOT EACH DATASET SEPARATELY
###################################

all_x <- c("minDistMidToBD_log10", "distMidToMid_log10", "minDistStartToBD_log10", "distStartToMid_log10", "minDistEndToBD_log10", "distEndToMid_log10")
y_var = "gene_logFC"

for(x_var in all_x) {
  foo <- foreach(ds = all_datasets) %dopar% {
    
    sub_data <- all_gene_tad_expr_dt[all_gene_tad_expr_dt$dataset == ds,]
  
    myx <- sub_data[,paste0(x_var)]
    myy <- sub_data[,paste0(y_var)]
    outFile <- file.path(outFolder, paste0(gsub(" - ", "_", ds), "_", y_var, "_vs_", x_var, "_densplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = myx,
      y = myy,
      cex.axis = axisCex,
      cex.lab = axisCex,
      xlab = paste0(gsub("_", " ", gsub("_log10", "", x_var), " [log10]")),
      ylab = paste0(gsub("_", " " , y_var)),
      main = paste0(y_var, " vs. ", x_var)
    )
    mtext(side=3, text = paste0(ds), font=3)
    mtext(side=3, text = paste0("(n=", nrow(sub_data),")"), font=3, line=-1)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }

  ###################################
  ################################### PLOT FOR ALL DATASETS
  ###################################
  myx <- all_gene_tad_expr_dt[,paste0(x_var)]
  myy <- all_gene_tad_expr_dt[,paste0(y_var)]
  outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x = myx,
    y = myy,
    cex.axis = axisCex,
    cex.lab = axisCex,
    xlab = paste0(gsub("_", " ", gsub("_log10", "", x_var)), " [log10]"),
    ylab = paste0(gsub("_", " " , y_var)),
    main = paste0(y_var, " vs. ", x_var)
  )
  mtext(side=3, text = paste0(nDS, " DS (n=", nrow(all_gene_tad_expr_dt),")"), font=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  

  myx <- all_gene_tad_expr_dt[,paste0(x_var)]
  myy <- all_gene_tad_expr_dt[,paste0(y_var)]
  outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_plotColType.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(
    x = myx,
    y = myy,
    cex.axis = axisCex,
    cex.lab = axisCex,
    xlab = paste0(gsub("_", " ", gsub("_log10", "", x_var)), " [log10]"),
    ylab = paste0(gsub("_", " " , y_var)),
    main = paste0(y_var, " vs. ", x_var),
    col = all_gene_tad_expr_dt$subtype_col,
    pch=16,
    cex=0.7
  )
  addSubtypeLeg(mypos="bottomleft", bty="n")
  mtext(side=3, text = paste0(nDS, " DS (n=", nrow(all_gene_tad_expr_dt),")"), font=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}

###################################
################################### BOXPLOT SIGNIF. AND DIST TO TAD
###################################
geneSignifThresh <- 0.05

all_gene_tad_expr_dt$geneSignif <- ifelse(all_gene_tad_expr_dt$gene_adjPval <= geneSignifThresh, "signif", "not.signif")

all_x <- c("minDistMidToBD_log10", "distMidToMid_log10", "minDistStartToBD_log10", "distStartToMid_log10", "minDistEndToBD_log10", "distEndToMid_log10")

curr_var = all_x[1]

for(curr_var in all_x) {
  
  outFile <- file.path(outFolder, paste0("allDS", "_", curr_var, "_", "geneSignif_", "boxplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  
  
  boxplot(formula(paste0(curr_var, "~", "geneSignif")), data = all_gene_tad_expr_dt, 
          ylab =curr_var,
          xlab="gene signif.")

  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  ###################################
  ################################### BOXPLOT SIGNIF. AND DIST TO TAD
  ###################################
  tad_agg_dt <- aggregate(formula(paste0(gsub("_log10", "", curr_var), "~ ", "region+hicds+exprds+TAD_adjCombPval")), 
                          data = all_gene_tad_expr_dt, FUN=mean)
  stopifnot(is.numeric(tad_agg_dt$TAD_adjCombPval))
  stopifnot(is.numeric(tad_agg_dt[ , paste0(gsub("_log10", "", curr_var))]))
  newColName <- paste0("TADmean_", gsub("_log10", "", curr_var))
  colnames(tad_agg_dt)[colnames(tad_agg_dt) == gsub("_log10", "", curr_var)] <- newColName
  
  
  outFile <- file.path(outFolder, paste0("allDS", "_meanTAD", curr_var, "_", "TADsignif_", "densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x=-log10(tad_agg_dt[,"TAD_adjCombPval"]),
    y=log10(tad_agg_dt[, newColName]),
    xlab = paste0("TAD_adjCombPval", " [-log10]"),
    ylab = paste0(newColName, " [log10]")
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}







##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

options(scipen=100)

# Rscript features_signif_TADs.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript features_signif_TADs.R   # to run all datasets in one shot

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "features_signif_TADs.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

buildData <- FALSE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
require(ggplot2)
require(reshape2)
require(metap)
require(lattice)
require(RColorBrewer)
hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("utils_fct.R")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script8c_name <- "8c_runAllDown"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script19_name <- "19_SAM_emp_measurement"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

samInFolder <- file.path("SAM_EMP_MEASUREMENT_MEANCORR")
stopifnot(dir.exists(samInFolder))

pipFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

corrPvalsFolder <- file.path("SAMPLE_MEANCORR_EMPPVALS")
stopifnot(dir.exists(corrPvalsFolder))

outFolder <- "FEATURES_SIGNIF_TADS"
dir.create(outFolder, recursive=TRUE)

corr_type <- "meanCorr"
samp_type <- "sameNbr"
permut_type <- "allDS"
fixSizeKb <- ""
sampName <- paste0("empPval-", samp_type, "-", corr_type, " - ", permut_type)

FDRthresh_seq <- seq(from=0.1, to=0.5, by=0.1)
pvalThresh_seq <- seq(from=0.01, to=0.05, by = 0.01)

myHeightGG <- length(pvalThresh_seq)*1.2
myWidthGG <- length(FDRthresh_seq)*1.2

twoSidedStouffer <- FALSE

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


if(buildData) {

  ### BUILD RATIO DOWN
  cat("... start building ratiodown and FC \n")
  allDS_rD_FC_DT <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds_ratioDown_DT <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      rd_file <- file.path(pipOutFolder, hicds, exprds, script8c_name, "all_obs_ratioDown.Rdata")
      stopifnot(file.exists(rd_file))
      tad_rD <- eval(parse(text = load(rd_file)))
      all_regs <- names(tad_rD)
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
      absMaxLogFC <- sapply(all_regs, function(reg) {
        g2t_tad_id <- g2t_DT$entrezID[g2t_DT$region==reg]
        de_tad_id <- names(geneList)[geneList %in% g2t_tad_id]
        stopifnot(de_tad_id %in% de_DT$genes)
        max(abs(de_DT$logFC[de_DT$genes %in% de_tad_id]))
      })
      names(absMaxLogFC) <- all_regs
      sdLogFC <- sapply(all_regs, function(reg) {
        g2t_tad_id <- g2t_DT$entrezID[g2t_DT$region==reg]
        de_tad_id <- names(geneList)[geneList %in% g2t_tad_id]
        stopifnot(de_tad_id %in% de_DT$genes)
        sd(de_DT$logFC[de_DT$genes %in% de_tad_id])
      })    
      names(sdLogFC) <- all_regs
      stopifnot(names(tad_rD) == names(sdLogFC))
      stopifnot(names(tad_rD) == names(absMaxLogFC))
      data.frame(
        hicds = hicds,
        exprds = exprds,
        region = names(tad_rD),
        ratioDown = as.numeric(tad_rD),
        maxAbsLogFC = as.numeric(absMaxLogFC),
        sdLogFC = as.numeric(sdLogFC),
        stringsAsFactors = FALSE
      )
    } # end-foreach iterating exprds
  } # end-foreach iterating hicds
  
  
  outFile <- file.path(outFolder, "allDS_rD_FC_DT.Rdata")  
  save(allDS_rD_FC_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  stopifnot(nrow(allDS_rD_FC_DT) > 0)
  
  ### BUILD SIGNIF ALONG FDR THRESH
  cat("... start building signif. along FDR thresh\n")
  allDS_signifFDR_DT <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds_signifFDR_DT <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      # PREPARE logFC and meanCorr observed data
      fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
      corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
      stopifnot(file.exists(fc_file))
      stopifnot(file.exists(corr_file))  
      tad_logFC <- eval(parse(text = load(fc_file)))
      tad_meanCorr <- eval(parse(text = load(corr_file)))
      all_regs <- names(tad_logFC)
      stopifnot(setequal(all_regs, names(tad_meanCorr)))
      tad_logFC <- tad_logFC[all_regs]
      tad_meanCorr <- tad_meanCorr[all_regs]
      # RETRIEVE FDR DATA FOR MEAN LOGFC
      logFC_FDR_file <- file.path(pipOutFolder, hicds, exprds, script19_name, "empFDR_list.Rdata")
      stopifnot(file.exists(logFC_FDR_file))
      all_FDR <- eval(parse(text = load(logFC_FDR_file)))
      logFC_FDR <- all_FDR[["empFDR_logFC"]]  # the names is the FC threshold, the value is the FDR
      stopifnot(length(logFC_FDR) > 0)
      # RETRIEVE FDR DATA FOR MEAN CORR
      meanCorr_FDR_file <-  file.path(samInFolder, samp_type, fixSizeKb, hicds, exprds, "all_empFDR.Rdata")
      stopifnot(file.exists(meanCorr_FDR_file))
      all_corr_FDR <- eval(parse(text = load(meanCorr_FDR_file)))
      meanCorr_FDR <- all_corr_FDR[[paste0("sample_", corr_type, "_", permut_type)]][["empFDR"]]  # the names is the meanCorr threshold, the value is the FDR
      stopifnot(length(meanCorr_FDR) > 0)
      # => SIGNIF TADs FOR VARIOUS FDR THRESHOLD
      # for each of the FDR threshold => FC cut-off, meanCorr cut-off => signif TADs
      cutoff_fdr <- FDRthresh_seq[1]
      signif_tads_dt <- foreach(cutoff_fdr = FDRthresh_seq , .combine='rbind') %do% {
        logFC_cut_off <- min(as.numeric(as.character(na.omit(names(logFC_FDR)[logFC_FDR <= cutoff_fdr]))))  # the smallest FC cut-off that leads to desired FDR; if not returns Inf
        meanCorr_cut_off <- min(as.numeric(as.character(na.omit(names(meanCorr_FDR)[meanCorr_FDR <= cutoff_fdr]))))
        stopifnot(names(tad_logFC) == names(tad_meanCorr))
        signif_tads <- (tad_logFC >= logFC_cut_off & tad_meanCorr >= meanCorr_cut_off)
        stopifnot(names(signif_tads) == names(tad_meanCorr))
        stopifnot(names(signif_tads) == all_regs)
        data.frame(
          hicds = hicds,
          exprds = exprds,
          region = names(signif_tads),
          thresh_FDR = cutoff_fdr,
          signif_FDR = as.logical(signif_tads), # to avoid having rownames
          stringsAsFactors = FALSE
        )
      } # end-foreach iterating over FDR threshs
    signif_tads_dt
    } # end-foreach iterating over exprds
    exprds_signifFDR_DT
  } # end-foreach iterating over hicds
  
  outFile <- file.path(outFolder, "allDS_signifFDR_DT.Rdata")  
  save(allDS_signifFDR_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  stopifnot(nrow(allDS_signifFDR_DT) > 0 )
  nIDs <- unique(paste0(allDS_signifFDR_DT$hicds, allDS_signifFDR_DT$exprds, allDS_signifFDR_DT$region))
  stopifnot(length(nIDs) == nrow(allDS_rD_FC_DT))
  
  ### BUILD SIGNIF ALONG PVAL THRESH
  cat("... start building signif. along pval thresh\n")
  allDS_signifPval_DT <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds_signifPval_DT <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      # PREPARE ADJ. COMBINED EMP PVAL      
      empPval_logFC_file <- file.path(pipOutFolder, hicds, exprds, script9_name, "emp_pval_meanLogFC.Rdata" )
      stopifnot(file.exists(empPval_logFC_file))
      empPval_logFC <- eval(parse(text = load(paste0(empPval_logFC_file))))
      all_regs <- names(empPval_logFC)
      empPval_meanCorr_file <- file.path(corrPvalsFolder, hicds, exprds, "all_cors_empPval_dt.Rdata")
      stopifnot(file.exists(empPval_meanCorr_file))
      empPval_meanCorr_dt <- eval(parse(text = load(paste0(empPval_meanCorr_file))))
      stopifnot(setequal(rownames(empPval_meanCorr_dt), all_regs))
      empPval_meanCorr_dt <- empPval_meanCorr_dt[all_regs,]
      stopifnot(sampName %in% colnames(empPval_meanCorr_dt))
      empPval_meanCorr <- setNames(empPval_meanCorr_dt[,paste0(sampName)], rownames(empPval_meanCorr_dt))
      stopifnot(names(empPval_meanCorr) == names(empPval_logFC))
      stopifnot(names(empPval_meanCorr) == all_regs)
      
      empPval_meanCorr <- empPval_meanCorr[all_regs]
      empPval_logFC <- empPval_logFC[all_regs]
      
      # COMPUTE THE COMBINED PVAL      
      empPval_comb <- sapply(all_regs, function(reg) {
        stopifnot(reg %in% names(empPval_logFC))
        stopifnot(reg %in% names(empPval_meanCorr))
        pval_logFC <- empPval_logFC[paste0(reg)]
        stopifnot(is.numeric(pval_logFC))
        pval_meanCorr <- empPval_meanCorr[paste0(reg)]
        stopifnot(is.numeric(pval_meanCorr))
        stouffer(c(pval_logFC, pval_meanCorr), two.tails=twoSidedStouffer)
      })
      stopifnot(names(empPval_comb) == all_regs)
      # ADJUST THE PVAL
      adj_empPval_comb <- p.adjust(empPval_comb, method="BH")
      stopifnot(names(adj_empPval_comb) == all_regs)
      # => SIGNIF TADs FOR VARIOUS PVAL THRESHOLD
      # for each of the threshold of empPvals -> retrieve signif TADs
      cutoff_pval=pvalThresh_seq[1]
      
      signif_tads_dt <- foreach(cutoff_pval = pvalThresh_seq , .combine='rbind') %do% {
        signif_tads_pval <- (adj_empPval_comb <= cutoff_pval)
        stopifnot(names(signif_tads_pval) == names(adj_empPval_comb))
        stopifnot(names(signif_tads_pval) == all_regs)
        data.frame(
          hicds = hicds,
          exprds = exprds,
          region = names(signif_tads_pval),
          thresh_pval = cutoff_pval,
          signif_pval = as.logical(signif_tads_pval), # to avoid having rownames
          stringsAsFactors = FALSE
        )
      } # end-foreach iterating over FDR threshs
      signif_tads_dt
    } # end-foreach iterating over exprds
    exprds_signifPval_DT
  } # end-foreach iterating over hicds
  
  
  outFile <- file.path(outFolder, "allDS_signifPval_DT.Rdata")  
  save(allDS_signifPval_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  stopifnot(nrow(allDS_signifPval_DT) > 0 )
  
  nIDs <- unique(paste0(allDS_signifPval_DT$hicds, allDS_signifPval_DT$exprds, allDS_signifPval_DT$region))
  stopifnot(length(nIDs) == nrow(allDS_rD_FC_DT))
  
  
  
  cat("... merge signifFDR <-> ratioDown\n")
  all_dt_signifFDR <- merge(allDS_signifFDR_DT, allDS_rD_FC_DT, by =c("hicds", "exprds", "region"), all = TRUE)
  stopifnot(nrow(all_dt_signifFDR) == nrow(allDS_signifFDR_DT))
  
  outFile <- file.path(outFolder, "all_dt_signifFDR.Rdata")  
  save(all_dt_signifFDR, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  cat("... merge signifFDR <-> ratioDown\n")
  all_dt_signifPval <- merge(allDS_signifPval_DT, allDS_rD_FC_DT, by =c("hicds", "exprds", "region"), all = TRUE)
  stopifnot(nrow(all_dt_signifPval) == nrow(allDS_signifPval_DT))
  
  outFile <- file.path(outFolder, "all_dt_signifPval.Rdata")  
  save(all_dt_signifPval, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
} else { # end-if build data
  outFile <- file.path(outFolder, "all_dt_signifPval.Rdata")  
  all_dt_signifPval <- eval(parse(text = load(outFile)))
  
  outFile <- file.path(outFolder, "all_dt_signifFDR.Rdata")  
  all_dt_signifFDR <- eval(parse(text = load(outFile)))
  
}

### PLOT THE SIGNIF PVAL FEATURES

all_dt_signifPval$adjCombPval <- ifelse(all_dt_signifPval$signif_pval, "signif.", "not signif.")
all_dt_signifPval$thresh_pval_label <- paste0("P=", all_dt_signifPval$thresh_pval)
all_vars <- c("ratioDown", "maxAbsLogFC", "sdLogFC")
plot_var=all_vars[1]
for(plot_var in all_vars) {
  outFile <- file.path(outFolder, paste0("allDS_signifPval_", plot_var, "_signif_vs_notsignif_density_lattice.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*3))
  myplot <- densityplot( formula(paste0("~",plot_var, "| thresh_pval_label")), groups = adjCombPval, data = all_dt_signifPval, #auto.key = TRUE, 
              # par.strip.text=list(cex=1), # width of the strip bar
              par.strip.text = list(cex = 1, font = 4, col = "brown"),
              layout = c(5, 5),
              scales=list(y=list(relation="free"),
                          x=list(relation="free")
              ),
              auto.key=list(title="", space = "bottom", cex=1.0, columns=length(unique(all_dt_signifPval$adjCombPval))),
              main = paste0(plot_var))
  print(myplot)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

### PLOT THE SIGNIF FDR

all_dt_signifFDR$FDR <- ifelse(all_dt_signifFDR$signif_FDR, "signif.", "not signif.")
all_dt_signifFDR$thresh_FDR_label <- paste0("FDR=", all_dt_signifPval$thresh_FDR)
all_vars <- c("ratioDown", "maxAbsLogFC", "sdLogFC")
plot_var=all_vars[1]
for(plot_var in all_vars) {
  outFile <- file.path(outFolder, paste0("allDS_signifFDR_", plot_var, "_signif_vs_notsignif_density_lattice.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*3))
  myplot <- densityplot( formula(paste0("~",plot_var, "| thresh_FDR_label")), groups = FDR, data = all_dt_signifFDR, #auto.key = TRUE, 
                         # par.strip.text=list(cex=1), # width of the strip bar
                         par.strip.text = list(cex = 1, font = 4, col = "brown"),
                         layout = c(5, 5),
                         scales=list(y=list(relation="free"),
                                     x=list(relation="free")
                         ),
                         auto.key=list(title="", space = "bottom", cex=1.0, columns=length(unique(all_dt_signifFDR$FDR))),
                         main = paste0(plot_var))
  print(myplot)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


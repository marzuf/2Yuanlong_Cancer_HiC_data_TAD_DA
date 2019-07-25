options(scipen=100)

# Rscript create_empPvalComb.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript create_empPvalComb.R   # to run all datasets in one shot

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "create_empPvalComb.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
require(ggplot2)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("utils_fct.R")

script0_name <- "0_prepGeneData"
script9_name <- "910000_runEmpPvalMeanTADLogFC"

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

outFolder <- "CREATE_EMPPVALCOMB"
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



### BUILD SIGNIF ALONG PVAL THRESH
cat("... start building signif. along pval thresh\n")
foo <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  foo <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
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
    meanCorr_meanLogFC_notAdjCombEmpPval <- sapply(all_regs, function(reg) {
      stopifnot(reg %in% names(empPval_logFC))
      stopifnot(reg %in% names(empPval_meanCorr))
      pval_logFC <- empPval_logFC[paste0(reg)]
      stopifnot(is.numeric(pval_logFC))
      pval_meanCorr <- empPval_meanCorr[paste0(reg)]
      stopifnot(is.numeric(pval_meanCorr))
      stouffer(c(pval_logFC, pval_meanCorr), two.tails=twoSidedStouffer)
    })
    stopifnot(names(meanCorr_meanLogFC_notAdjCombEmpPval) == all_regs)

    outFile <- file.path(outFolder, samp_type, fixSizeKb, hicds, exprds, paste0(corr_type, "_", "meanLogFC", "_notAdjCombEmpPval.Rdata"))
    dir.create(dirname(outFile), recursive = TRUE)
    save(meanCorr_meanLogFC_notAdjCombEmpPval, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    
  } # end-foreach iterating over exprds
} # end-foreach iterating over hicds


  



##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


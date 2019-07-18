options(scipen=100)

# Rscript combined_pvals_meth_cmp.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript combined_pvals_meth_cmp.R   # to run all datasets in one shot

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "combined_pvals_meth_cmp.R"

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

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("utils_fct.R")

script0_name <- "0_prepGeneData"
script9_name <- "9_runEmpPvalMeanTADLogFC"
# ../Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/9_runEmpPvalMeanTADLogFC/emp_pval_meanLogFC.Rdata

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4


pipOutFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

corrPvalsFolder <- file.path("SAMPLE_MEANCORR_EMPPVALS")
stopifnot(dir.exists(corrPvalsFolder))

outFolder <- "COMBINED_PVALS_METH_CMP"
dir.create(outFolder, recursive=TRUE)

corr_type <- "meanCorr"
samp_type <- "sameNbr"
permut_type <- "allDS"

sampName <- paste0("empPval-", samp_type, "-", corr_type, " - ", permut_type)

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
  allDS_combPvals_DT <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    
    exprds_combPvals_DT <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start ", hicds, " - " , exprds, "\n")
  
      
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
      
      stouffer2tails_pvals <- sapply(all_regs, function(reg) {
        stopifnot(reg %in% names(empPval_logFC))
        stopifnot(reg %in% names(empPval_meanCorr))
        pval_logFC <- empPval_logFC[paste0(reg)]
        stopifnot(is.numeric(pval_logFC))
        pval_meanCorr <- empPval_meanCorr[paste0(reg)]
        stopifnot(is.numeric(pval_meanCorr))
        s1 <- stouffer(c(pval_logFC, pval_meanCorr), two.tails=TRUE)
        s2 <- numeric(sumz(c(pval_logFC, pval_meanCorr)/2)$p) * 2
        stopifnot(abs(s1-s2) < 10^-5)
        s1
      })
      stopifnot(names(stouffer2tails_pvals) == all_regs)
      
      
      all_combined_pvals <- lapply(all_regs, function(reg) {
        stopifnot(reg %in% names(empPval_logFC))
        stopifnot(reg %in% names(empPval_meanCorr))
        pval_logFC <- empPval_logFC[paste0(reg)]
        stopifnot(is.numeric(pval_logFC))
        pval_meanCorr <- empPval_meanCorr[paste0(reg)]
        stopifnot(is.numeric(pval_meanCorr))
        all_pvals <- allmetap(c(pval_meanCorr, pval_logFC), method="all")$p
        setNames(as.numeric(all_pvals), names(all_pvals))
      })
      names(all_combined_pvals) <- all_regs
      
      all_combPvals_dt <- as.data.frame(do.call(rbind, all_combined_pvals))
      
      all_combPvals_dt$sumzTwoSided <- stouffer2tails_pvals
      
      all_combPvals_dt$hicds <- hicds
      all_combPvals_dt$exprds <- exprds
      all_combPvals_dt$region <- rownames(all_combPvals_dt)
      rownames(all_combPvals_dt) <- NULL
      all_combPvals_dt
    } # end-foreach iterating over exprds
    exprds_combPvals_DT
  } # end-foreach iterating over hicds  
    
  outFile <- file.path(outFolder, "allDS_combPvals_DT.Rdata")  
  save(allDS_combPvals_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "allDS_combPvals_DT.Rdata")  
  allDS_combPvals_DT <- eval(parse(text = load(outFile)))
}
# allDS_combPvals_DT

nDS <- length(unique(paste0(allDS_combPvals_DT$hicds,"-", allDS_combPvals_DT$exprds)))

x_var <- "sumz"
stopifnot(x_var %in% colnames(allDS_combPvals_DT))
my_x <- allDS_combPvals_DT[,paste0(x_var)]
all_y <- colnames(allDS_combPvals_DT)[! colnames(allDS_combPvals_DT) %in% c(x_var, "region", "exprds", "hicds")]
y_var=all_y[1]
for(y_var in all_y) {
  stopifnot(y_var %in% colnames(allDS_combPvals_DT))
  my_y <- allDS_combPvals_DT[,paste0(y_var)]

  if(all(is.na(my_y))) next
      
  outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_pvalComb_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myHeight))
  densplot(
    x = my_x,
    xlab = paste0(x_var),
    y = my_y,
    ylab = paste0(y_var),
    main = paste0("Comb. p-val.: ", y_var, " vs. ", x_var),
    cex = 0.7,
    cex.lab = axisCex,
    cex.axis = axisCex
  )
  mtext(side = 3, text = paste0("nDS = ", nDS, ";  one-sided"), font = 3)
  addCorr(x=my_x, legPos="topleft",
          y=my_y, bty='n')
  curve(1*x, add=TRUE, lty=2, col="darkgrey")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_pvalComb_densplot_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myHeight))
  densplot(
    x = -log10(my_x),
    xlab = paste0(x_var, " [-log10]"),
    y = -log10(my_y),
    ylab = paste0(y_var, " [-log10]"),
    main = paste0("Comb. p-val.: ", y_var, " vs. ", x_var),
    cex = 0.7,
    cex.lab = axisCex,
    cex.axis = axisCex
  )
  mtext(side = 3, text = paste0("nDS = ", nDS, ";  one-sided"), font = 3)
  addCorr(x=my_x, legPos="topleft",
          y=my_y, bty='n')
  curve(1*x, add=TRUE, lty=2, col="darkgrey")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}# end-iterating over y_var
    


##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


# [...] we would want to keep track of the direction of the effect when we combine the information from the different studies. 
# The solution is to calculate the P-value from each test as a one-sided P-value, with P-values close to zero having a consistent meaning and P-values close to one having the opposite meaning. 
# After combining the P-values, if desired the resulting combined P can be again *converted to a two-tailed test by multiplying it by two.*

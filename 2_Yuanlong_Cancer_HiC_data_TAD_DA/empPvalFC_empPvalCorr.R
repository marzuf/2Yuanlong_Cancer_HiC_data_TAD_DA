# Rscript empPvalFC_empPvalCorr.R

script_name <- "empPvalFC_empPvalCorr.R"
cat("> Start ", script_name, "\n")
startTime <- Sys.time()

require(foreach)
require(doMC)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

registerDoMC(40)

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight

cexPlot <- 1.2


outFold <- "EMPPVALFC_EMPPVALCORR"
dir.create(outFold)

dataFolder <- "COEXPR_BETWEEN_WITHIN_ALL"

pipFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

corrAcrossFolder <- "SAMPLE_MEANCORR_EMPPVALS"
stopifnot(dir.exists(corrAcrossFolder))

# all_pvalComb_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_combined.Rdata", full.names = FALSE)
# stopifnot(length(all_pvalComb_files) > 0)

all_fc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_meanLogFC.Rdata", full.names = FALSE)
stopifnot(length(all_fc_files) > 0)

all_meanCorr_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_meanCorr.Rdata", full.names = FALSE)
stopifnot(length(all_meanCorr_files) > 0)
# 
# dataFile <- file.path(dataFolder, "allData_within_between_coexpr.Rdata")
# stopifnot(file.exists(dataFile))
# allData_within_between_coexpr <- eval(parse(text = load(dataFile)))

# all_ratioDown_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_obs_ratioDown.Rdata", full.names = FALSE)
# stopifnot(length(all_ratioDown_files) > 0)

all_meanCorrAcross_files <-  list.files(corrAcrossFolder  , recursive = TRUE, pattern="all_cors_empPval_dt.Rdata", full.names = FALSE)
stopifnot(length(all_meanCorrAcross_files) > 0)

### BUILD THE LOGFC TABLE
cat("... start build fc_DT \n")
fc_file = all_fc_files[1]
fc_DT <- foreach(fc_file = all_fc_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, fc_file)
  stopifnot(file.exists(curr_file))
  tad_fc <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(fc_file))
  hicds <- dirname(dataset)
  exprds <- basename(dataset)
  data.frame(
    dataset = dataset,
    hicds = hicds,
    exprds=exprds,
    region = names(tad_fc),
    empPval_meanFC = as.numeric(tad_fc),
    stringsAsFactors = FALSE
  )
}
### BUILD THE MEANCORR TABLE
cat("... start build meanCorr_DT \n")
meanCorr_file = all_meanCorr_files[1]
meanCorr_DT <- foreach(meanCorr_file = all_meanCorr_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, meanCorr_file)
  stopifnot(file.exists(curr_file))
  tad_meanCorr <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(meanCorr_file))
  hicds <- dirname(dataset)
  exprds <- basename(dataset)
  data.frame(
    dataset = dataset,
    hicds = hicds,
    exprds=exprds,
    region = names(tad_meanCorr),
    empPval_meanCorr = as.numeric(tad_meanCorr),
    stringsAsFactors = FALSE
  )
}


### BUILD THE MEANCORR TABLE
cat("... start build meanCorrAcross_DT \n")
meanCorrAcross_file = all_meanCorrAcross_files[1]
meanCorrAcross_DT <- foreach(meanCorrAcross_file = all_meanCorrAcross_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(corrAcrossFolder, meanCorrAcross_file)
  stopifnot(file.exists(curr_file))
  meanCorrDT <- eval(parse(text = load(curr_file)))
  dataset <- dirname(meanCorrAcross_file)
  hicds <- dirname(dataset)
  exprds <- basename(dataset)
  
  meanCorrDT$region <-rownames(meanCorrDT)
  rownames(meanCorrDT) <- NULL
  meanCorrDT$dataset <- dataset
  meanCorrDT$hicds <- hicds
  meanCorrDT$exprds <- exprds
  
  meanCorrDT
}

stopifnot(nrow(meanCorr_DT) == nrow(fc_DT))
stopifnot(nrow(meanCorr_DT) == nrow(meanCorrAcross_DT))


id_cols <-  c("dataset", "hicds", "exprds", "region")

all_dt <- merge(meanCorr_DT, fc_DT, by = id_cols, all = TRUE)
stopifnot(nrow(all_dt) == nrow(fc_DT))


all_dt <- merge(all_dt, meanCorrAcross_DT, by =  id_cols, all = TRUE)
stopifnot(nrow(all_dt) == nrow(meanCorrAcross_DT))


outFile <- file.path(outFold, "all_dt.Rdata")
save(all_dt, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


nDS <- length(unique(all_dt$dataset))
nTADs <- nrow(all_dt)

#############################################################################################################################
############################################################################################################################# ... vs. empPvalMeanFC & ... vs. empPvalMeanCorr
#############################################################################################################################

all_x <- c("empPval_meanCorr", "empPval_meanFC")

all_y <- colnames(all_dt)[!colnames(all_dt) %in% c(all_x, id_cols)]

for(x_var in all_x) {
  
  stopifnot(x_var %in% colnames(all_dt))
  
  myx <- all_dt[,paste0(x_var)]
  
  foo <- foreach(y_var = all_y) %dopar% {
    
    
    
    outFile <- file.path(outFold, paste0(y_var, "_vs_", x_var, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = myx,
      y = all_dt[,paste0(y_var)],
      cex.axis = cexPlot,
      cex.lab = cexPlot,
      xlab = paste0(x_var),
      ylab = paste0(y_var),
      main = paste0(y_var, " vs.\n", x_var)
    )
    mtext(side = 3, text = paste0("nDS = ", nDS, " - nTADs = ", nTADs), font = 3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    outFile <- file.path(outFold, paste0(y_var, "_vs_", x_var, "_log10.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = - log10(myx),
      y = - log10(all_dt[,paste0(y_var)]),
      cex.axis = cexPlot,
      cex.lab = cexPlot,
      xlab = paste0(x_var, " [-log10]"),
      ylab = paste0(y_var, " [-log10]"),
      main = paste0(y_var, " vs.\n", x_var)
    )
    mtext(side = 3, text = paste0("nDS = ", nDS, " - nTADs = ", nTADs), font = 3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
  
  
}


#############################################################################################################################
############################################################################################################################# empPvalMeanCorr vs. empPvalMeanFC
#############################################################################################################################

all_x <- c("empPval_meanFC")

all_y <- c("empPval_meanCorr")

for(x_var in all_x) {
  
  stopifnot(x_var %in% colnames(all_dt))
  
  myx <- all_dt[,paste0(x_var)]
  
  foo <- foreach(y_var = all_y) %dopar% {
    
    
    
    outFile <- file.path(outFold, paste0(y_var, "_vs_", x_var, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = myx,
      y = all_dt[,paste0(y_var)],
      cex.axis = cexPlot,
      cex.lab = cexPlot,
      xlab = paste0(x_var),
      ylab = paste0(y_var),
      main = paste0(y_var, " vs.\n", x_var)
    )
    mtext(side = 3, text = paste0("nDS = ", nDS, " - nTADs = ", nTADs), font = 3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    outFile <- file.path(outFold, paste0(y_var, "_vs_", x_var, "_log10.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = - log10(myx),
      y = - log10(all_dt[,paste0(y_var)]),
      cex.axis = cexPlot,
      cex.lab = cexPlot,
      xlab = paste0(x_var, " [-log10]"),
      ylab = paste0(y_var, " [-log10]"),
      main = paste0(y_var, " vs.\n", x_var)
    )
    mtext(side = 3, text = paste0("nDS = ", nDS, " - nTADs = ", nTADs), font = 3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
  
  
}

#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))

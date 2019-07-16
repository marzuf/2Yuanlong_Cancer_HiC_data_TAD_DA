# as in the Tusher article, plott the observed FC vs. the average FC from permut

options(scipen=100)

# Rscript scatterplot_FDR_byDS.R

script_name <- "scatterplot_FDR_byDS.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(tools)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4


outFolder <- "SCATTERPLOT_FDR_BYDS"
dir.create(outFolder, recursive=TRUE)

pipOutFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

all_fc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanLogFC_TAD.Rdata", full.names = FALSE)
stopifnot(length(all_fc_files) > 0)

all_fcPermut_files <- list.files(pipOutFolder, recursive = TRUE, pattern="meanLogFC_permDT.Rdata", full.names = FALSE)
stopifnot(length(all_fcPermut_files) > 0)

script6_name <- "6_runPermutationsMeanLogFC"
# /6_runPermutationsMeanLogFC/meanLogFC_permDT.Rdata"

stopifnot(length(all_fc_files) == length(all_fcPermut_files) )


### BUILD THE LOGFC TABLE
fc_file = all_fc_files[1]
fc_DT <- foreach(fc_file = all_fc_files, .combine = 'rbind') %dopar% {
  dataset <- dirname(dirname(fc_file))  
  cat("... ", dataset, " - build obs. meanFC table \n")
  
  curr_file <- file.path(pipOutFolder, fc_file)
  stopifnot(file.exists(curr_file))
  tad_fc <- eval(parse(text = load(curr_file)))
  
  fc_DT <- data.frame(
    dataset = dataset,
    region = names(tad_fc),
    meanFC = as.numeric(tad_fc),
    stringsAsFactors = FALSE
  )
  
  permut_file <- file.path(pipOutFolder, dataset, script6_name, "meanLogFC_permDT.Rdata")
  stopifnot(file.exists(permut_file))
  tad_permDT <- eval(parse(text = load(permut_file)))
  tad_avgFC <- apply(tad_permDT, 1, mean)
  stopifnot(!is.na(tad_avgFC))
  
  cat("... ", dataset, " - build exp. meanFC table \n")
  
  permut_fc_DT <- data.frame(
    dataset = dataset,
    region = names(tad_avgFC),
    meanFC_avgPermut = as.numeric(tad_avgFC),
    stringsAsFactors = FALSE
  )
  
   stopifnot(setequal(permut_fc_DT$region, fc_DT$region))
   stopifnot(nrow(permut_fc_DT$region) == nrow(fc_DT$region))
    
   all_FC_DT <- merge(fc_DT, permut_fc_DT, by = c("dataset", "region"), all = TRUE)
   
   stopifnot(nrow(fc_DT) == nrow(all_FC_DT))
   stopifnot(nrow(permut_fc_DT) == nrow(all_FC_DT))
   
   ############################################################## meanFC densplot
   cat("... ", dataset, " - start plotting meanFC  \n")
   myy <- all_FC_DT$meanFC
   myx <- all_FC_DT$meanFC_avgPermut

  nTADs <- nrow(all_FC_DT)
   
   outFile <- file.path(outFolder,dataset, paste0("observed_expected_meanFC_densplot.", plotType))
   dir.create(dirname(outFile), recursive = TRUE)
   do.call(plotType, list(outFile, height=myHeight, width=myWidth))
   densplot(
     x = myx,
     y = myy,
     xlab = paste0("avg. permut. TAD meanFC"),
     ylab = paste0("obs. TAD meanFC"),
     cex.lab = plotCex,
     cex.axis = plotCex,
     main = paste0("obs. vs. expected TAD meanFC")
   )
   mtext(side = 3, text = paste0(dataset, "; nTADs = ", nTADs), font=3)
   curve(1*x, add=TRUE, lty=2, col="black", lwd=1.5)
   addCorr(
     x = myx, y = myy,
     legPos = "topleft", bty="n"
   )
   foo <- dev.off()
   cat(paste0("... written: ", outFile, "\n"))
   
   
   ############################################################## meanFC scatterplot grey threshold
   fc_threshold <- 1.5
   dotcols <- ifelse( abs(all_FC_DT$meanFC) >= fc_threshold, "black", "grey" )
   
   outFile <- file.path(outFolder, dataset, paste0("observed_expected_meanFC_threshold", fc_threshold, "grey_scatterplot.", plotType))
   dir.create(dirname(outFile), recursive = TRUE)
   do.call(plotType, list(outFile, height=myHeight, width=myWidth))
   plot(
     x = myx,
     y = myy,
     col = dotcols,
     pch = 16,
     cex = 0.7,
     xlab = paste0("avg. permut. TAD meanFC"),
     ylab = paste0("obs. TAD meanFC"),
     cex.lab = plotCex,
     cex.axis = plotCex,
     main = paste0("obs. vs. expected TAD meanFC")
   )
   mtext(side = 3, text = paste0(dataset, "; nTADs = ", nTADs), font=3)
   curve(1*x, add=TRUE, lty=2, col="black", lwd=1.5)
   addCorr(
     x = myx, y = myy,
     legPos = "topleft", bty="n"
   )
   foo <- dev.off()
   cat(paste0("... written: ", outFile, "\n"))
   
 
  
}





















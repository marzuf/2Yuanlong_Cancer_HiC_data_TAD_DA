# as in the Tusher article, plott the observed FC vs. the average FC from permut

options(scipen=100)

# Rscript scatterplot_FDR_allDS.R

script_name <- "scatterplot_FDR_allDS.R"

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

script6_name <- "6_runPermutationsMeanLogFC"

outFolder <- "SCATTERPLOT_FDR_ALLDS"
dir.create(outFolder, recursive=TRUE)

pipOutFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

all_fc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanLogFC_TAD.Rdata", full.names = FALSE)
stopifnot(length(all_fc_files) > 0)

all_fcPermut_files <- list.files(pipOutFolder, recursive = TRUE, pattern="meanLogFC_permDT.Rdata", full.names = FALSE)
stopifnot(length(all_fcPermut_files) > 0)

all_fcPermut_files <- all_fcPermut_files[grepl(script6_name, all_fcPermut_files)]

stopifnot(length(all_fc_files) == length(all_fcPermut_files) )


### BUILD THE LOGFC TABLE
cat("... build obs. meanFC table \n")
fc_file = all_fc_files[1]
fc_DT <- foreach(fc_file = all_fc_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, fc_file)
  stopifnot(file.exists(curr_file))
  tad_fc <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(fc_file))
  data.frame(
    dataset = dataset,
    region = names(tad_fc),
    meanFC = as.numeric(tad_fc),
    stringsAsFactors = FALSE
  )
}


### BUILD THE PERMUT AVERAGE LOGFC TABLE
cat("... build exp. meanFC table \n")
fc_file = all_fcPermut_files[1]
avgFcPermut_DT <- foreach(fc_file = all_fcPermut_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, fc_file)
  stopifnot(file.exists(curr_file))
  tad_permDT <- eval(parse(text = load(curr_file)))
  tad_avgFC <- apply(tad_permDT, 1, mean)
  stopifnot(!is.na(tad_avgFC))
  dataset <- dirname(dirname(fc_file))
  data.frame(
    dataset = dataset,
    region = names(tad_avgFC),
    meanFC_avgPermut = as.numeric(tad_avgFC),
    stringsAsFactors = FALSE
  )
}

stopifnot(nrow(fc_DT) == nrow(avgFcPermut_DT))

all_FC_DT <- merge(fc_DT, avgFcPermut_DT, by = c("dataset", "region"), all = TRUE)
stopifnot(nrow(fc_DT) == nrow(all_FC_DT))
stopifnot(nrow(avgFcPermut_DT) == nrow(all_FC_DT))

outFile <- file.path(outFolder, "all_FC_DT.Rdata")
save(all_FC_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

cat("... start plotting meanFC  \n")

############################################################## meanFC densplot
myy <- all_FC_DT$meanFC
myx <- all_FC_DT$meanFC_avgPermut

nDS <- length(unique(all_FC_DT$dataset))
nTADs <- nrow(all_FC_DT)

outFile <- file.path(outFolder, paste0("observed_expected_meanFC_densplot.", plotType))
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
mtext(side = 3, text = paste0("nDS = ", nDS, "; nTADs = ", nTADs), font=3)
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

outFile <- file.path(outFolder, paste0("observed_expected_meanFC_threshold", fc_threshold, "grey_scatterplot.", plotType))
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
mtext(side = 3, text = paste0("nDS = ", nDS, "; nTADs = ", nTADs), font=3)
curve(1*x, add=TRUE, lty=2, col="black", lwd=1.5)
addCorr(
  x = myx, y = myy,
  legPos = "topleft", bty="n"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))





















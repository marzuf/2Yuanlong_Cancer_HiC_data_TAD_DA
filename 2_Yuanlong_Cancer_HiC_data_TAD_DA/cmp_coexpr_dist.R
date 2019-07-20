# Rscript cmp_coexpr_dist.R

startTime <- Sys.time()

options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

script_name <- "cmp_coexpr_dist.R"

outFolder <- "CMP_COEXPR_DIST"
dir.create(outFolder)

all_sampTypes <- c("sameKb", "fixKb", "sameNbr")


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight * 1.2

nSampled <- 1000
alternativeCorr_file <- file.path("COEXPR_ALTERNATIVE_SAMPLINGS", nSampled, "ds_all_corr_data.Rdata")
stopifnot(file.exists(alternativeCorr_file))

samp_type <- "sameNbr"
fixKbSize <- ""
meanCorr_file <- file.path("COEXPR_AROUND_TADS", samp_type, 
                           fixKbSize, "all_ds_around_TADs_corr.Rdata")
stopifnot(file.exists(meanCorr_file))

pipFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
all_obsCorr_files <- list.files(pipOutFolder, pattern="all_meanCorr_TAD.Rdata", recursive=TRUE, full.names=TRUE)


### PREPARE ALTERNATIVE SAMPLING CORR. DATA
alternativeCorr_data <- eval(parse(text = load(alternativeCorr_file)))

alternativeCorr_values <- unlist(lapply(alternativeCorr_data, 
                                        function(sub1) lapply(sub1, function(sub2) lapply(sub2, function(x) x[["meanCorr_all"]]))))
alternativeCorr_DT <- data.frame(
  meanCorr = as.numeric(alternativeCorr_values),
  meanCorr_type = paste0("sampleV2_",nSampled),
  stringsAsFactors = FALSE
)

### PREPARE MEAN CORR DATA
sampCorr_data <- eval(parse(text = load(meanCorr_file)))
sampCorr_values <- unlist(lapply(sampCorr_data, function(sub_data) lapply(sub_data, function(x)x[["meanCorr"]])))
sampCorr_DT <- data.frame(
  meanCorr=as.numeric(sampCorr_values),
  meanCorr_type="sample_sameNbr",
  stringsAsFactors = FALSE
)

### PREPARE OBSERVED DATA
meanCorr_DT <- foreach(c_file = all_obsCorr_files, .combine='rbind') %dopar% {
  stopifnot(file.exists(c_file))
  tad_corr <- eval(parse(text = load(c_file)))
  data.frame(
#    region = names(tad_corr),
    meanCorr = tad_corr,
    meanCorr_type = "observed",
    stringsAsFactors = FALSE
  )
}


all_corr_DT <- rbind(alternativeCorr_DT, rbind(meanCorr_DT, sampCorr_DT))

outFile <- file.path(outFolder, paste0("all_corr_DT.Rdata"))
save(all_corr_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("meanCorr_cmp_multidens.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  split(all_corr_DT$meanCorr, all_corr_DT$meanCorr_type),
  plotTit = paste0("meanCorr comparison"),
  my_xlab = paste0("meanCorr")
  )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





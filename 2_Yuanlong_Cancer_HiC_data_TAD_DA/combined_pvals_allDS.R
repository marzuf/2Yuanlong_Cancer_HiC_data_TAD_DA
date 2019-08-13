options(scipen=100)

# Rscript combined_pvals_allDS.R 

setDir=""

script_name <- "combined_pvals_allDS.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
require(ggplot2)
require(reshape2)


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("utils_fct.R")

script0_name <- "0_prepGeneData"
script9_name <- "9_runEmpPvalMeanTADLogFC"
# ../Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/9_runEmpPvalMeanTADLogFC/emp_pval_meanLogFC.Rdata

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4

myHeightGG <- 7
myWidthGG <- myHeightGG 

pval_thresh <- 0.05

pipOutFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))


corrPvalsFolder <- file.path("SAMPLE_MEANCORR_EMPPVALS")
stopifnot(dir.exists(corrPvalsFolder))

adjcombPvalsFolder <- file.path("COMBINED_PVALS")
stopifnot(dir.exists(adjcombPvalsFolder))

outFolder <- "COMBINED_PVALS_ALLDS"
dir.create(outFolder, recursive=TRUE)

all_corrPval_files <- list.files(corrPvalsFolder, pattern = "all_cors_empPval_dt.Rdata", recursive = TRUE, full.names = TRUE)
stopifnot(length(all_corrPval_files) > 0)

all_adjCombPval_files <- list.files(adjcombPvalsFolder, pattern="all_adjCombinedPvals_dt.Rdata", recursive = TRUE, full.names = TRUE)
stopifnot(length(all_adjCombPval_files) == length(all_corrPval_files))

# COMBINE THE CORR PVAL VALUES
curr_file = all_corrPval_files[1]

all_corrPvals_dt <- foreach(curr_file = all_corrPval_files, .combine='rbind') %dopar% {
  stopifnot(file.exists(curr_file))    
  hicds <- basename(dirname(dirname(curr_file)))
  exprds <- basename(dirname(curr_file))
  corr_dt <- eval(parse(text = load(curr_file)))
  rownames(corr_dt) <- NULL
  corr_dt$dataset <- paste0(hicds, "-", exprds)
  corr_dt
}


# COMBINE THE ADJ. COMB PVAL VALUES
curr_file = all_adjCombPval_files[1]
all_adjCombPvals_dt <-  foreach(curr_file = all_adjCombPval_files, .combine='rbind') %dopar% {
  stopifnot(file.exists(curr_file))    
  hicds <- basename(dirname(dirname(curr_file)))
  exprds <- basename(dirname(curr_file))
  corr_dt <- eval(parse(text = load(curr_file)))
  corr_dt <- as.data.frame(corr_dt)
  rownames(corr_dt) <- NULL
  corr_dt$dataset <- paste0(hicds, "-", exprds)
  corr_dt
}


nDS <- length(unique(all_adjCombPvals_dt$dataset))
nDS2 <- length(unique(all_corrPvals_dt$dataset))

stopifnot(nDS == nDS2)





###################################################################################
################################################################################### DENSITY CORR PVALS
###################################################################################
all_corrPvals_dt$dataset <- NULL
corrPvals_list1 <- apply(all_corrPvals_dt, 2, as.list)
corrPvals_list <- lapply(corrPvals_list1, as.numeric)
corrPvals_list_log10 <- lapply(corrPvals_list, function(x) -log10(as.numeric(x)))

outFile <- file.path(outFolder, paste0("allCorr_corrEmpPval_densplot.", plotType))
dir.create(dirname(outFile), recursive = TRUE)
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
plot_multiDens(
  corrPvals_list_log10,
  plotTit = paste0("emp. p-val meanCorr - all DS "),
  my_xlab = paste0("emp. p-val [-log10]")
)
mtext(side=3, text=paste0("nDS = ", nDS), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

toKeep1 <- which(grepl("meanCorr - allDS", names(corrPvals_list_log10)))
toKeep2 <- which(grepl("meanCorrLeftRight - allDS", names(corrPvals_list_log10)))

outFile <- file.path(outFolder, paste0("meanCorr_meanCorrLeftRight_allDS_corrEmpPval_densplot.", plotType))
dir.create(dirname(outFile), recursive = TRUE)
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
plot_multiDens(
  corrPvals_list_log10[c(toKeep1, toKeep2)],
  plotTit = paste0("emp. p-val meanCorr - all DS "),
  my_xlab = paste0("emp. p-val [-log10]")
)
mtext(side=3, text=paste0("nDS = ", nDS), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

###################################################################################
################################################################################### DENSITY COMBINED PVALS
###################################################################################
all_adjCombPvals_dt_0 <- all_adjCombPvals_dt
all_adjCombPvals_dt$dataset <- NULL
combinedPvals_list1 <- apply(all_adjCombPvals_dt, 2, as.list)
combinedPvals_list <- lapply(combinedPvals_list1, as.numeric)
combinedPvals_list_log10 <- lapply(combinedPvals_list1, function(x) -log10(as.numeric(x)))

outFile <- file.path(outFolder,  paste0("allCorr_adjCombEmpPval_densplot.", plotType))
dir.create(dirname(outFile), recursive = TRUE)
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
plot_multiDens(
  combinedPvals_list_log10,
  plotTit = paste0("adj. emp. p-val combined - all DS "),
  my_xlab = paste0("emp. p-val [-log10]")
)
mtext(side=3, text=paste0("nDS = ", nDS), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


toKeep1 <- which(grepl("meanCorr - allDS", names(combinedPvals_list_log10)))
toKeep2 <- which(grepl("meanCorrLeftRight - allDS", names(combinedPvals_list_log10)))

outFile <- file.path(outFolder,paste0("meanCorr_meanCorrLeftRight_allDS_adjCombEmpPval_densplot.", plotType))
dir.create(dirname(outFile), recursive = TRUE)
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
plot_multiDens(
  combinedPvals_list_log10[c(toKeep1, toKeep2)],
  plotTit = paste0("adj. emp. p-val combined - all DS "),
  my_xlab = paste0("emp. p-val [-log10]")
)
mtext(side=3, text=paste0("nDS = ", nDS), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


###################################################################################
################################################################################### BOXPLOT NBR SIGNIF
###################################################################################
melt_comb_dt <- melt(all_adjCombPvals_dt_0, id ="dataset")
colnames(melt_comb_dt) <- c("dataset", "pvalType", "pval")
stopifnot(is.numeric(melt_comb_dt$pval))


plot_dt <- do.call(rbind, by(melt_comb_dt, melt_comb_dt$dataset, function(ds_dt) {
  do.call(rbind, by(ds_dt, ds_dt$pvalType, function(sub_dt) {
    ds <- unique(sub_dt$dataset)
    stopifnot(length(ds) == 1)
    corrType <- unique(sub_dt$pvalType)
    stopifnot(length(corrType) == 1)
    nSignif <- sum(sub_dt$pval <= pval_thresh)
    stopifnot(is.numeric(nSignif))
    data.frame(
      dataset = ds,
      corrType = corrType,
      nSignif = nSignif,
      stringsAsFactors=FALSE
    )
  }))
}))
plot_dt <- data.frame(plot_dt)
rownames(plot_dt) <- NULL

nSignif_dt <- plot_dt

nSignif_dt$samp_type <- gsub("adjComb_empPval-(.+)-mean.+", "\\1", nSignif_dt$corrType)
nSignif_dt$corr_type <- gsub("adjComb_empPval-.+-(mean.+)", "\\1", nSignif_dt$corrType)
nSignif_dt$corr_type <- gsub(" - ", "\n", nSignif_dt$corr_type)

p_nsignif <-  ggplot(nSignif_dt, aes(x = samp_type, y = nSignif, fill = samp_type)) + 
  geom_boxplot()+
  facet_grid(~corr_type, switch="x") + 
  coord_cartesian(expand = FALSE) +
  ggtitle(paste0("all DS (n = ", nDS, ")"), subtitle = paste0("# signif. TADs - adj. comb. p-val. thresh = ", pval_thresh) )+
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0("# signif. TADs"),
                     breaks = scales::pretty_breaks(n = 20))+
  # scale_fill_brewer(palette="YlOrRd")+
  labs(fill  = "") +
  theme( # Increase size of axis lines
    strip.text = element_text(size = 10),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    strip.text.x = element_text(size = 10),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .3, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )
outFile <- file.path(outFolder, paste0("nSignif_by_corrType_sampType_signifThresh", pval_thresh , "_barplot.", plotType))
dir.create(dirname(outFile), recursive = TRUE)
ggsave(p_nsignif, filename = outFile, height = myHeightGG, width = myWidthGG*1.4)
cat(paste0("... written: ", outFile, "\n"))








##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))







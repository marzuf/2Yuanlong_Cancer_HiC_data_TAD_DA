options(scipen=100)

script_name <- "look_sampling_features.R"
cat("... start ", script_name, "\n")
startTime <- Sys.time()


require(reshape2)
require(foreach)
require(doMC)
require(ggplot2)
require(hrbrthemes)
require(ggthemes)

registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

fixKbSize <- 1000000


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
myHeightGG <- 7
myWidthGG <- myHeightGG 

outFold <- "LOOK_SAMPLING_FEATURES"
dir.create(outFold, recursive = TRUE)

all_ds_sample_around_TADs_sameKb_file <- file.path("CREATE_SAMPLE_AROUND_TADS_SAMEKB_TWOSIDED", "all_ds_sample_around_TADs_sameKb.Rdata")
all_ds_sample_around_TADs_fixKb_file <- file.path("CREATE_SAMPLE_AROUND_TADS_FIXKB_TWOSIDED", fixKbSize, "all_ds_sample_around_TADs_fixKb.Rdata")
all_ds_sample_around_TADs_sameNbr_file <- file.path("CREATE_SAMPLE_AROUND_TADS_SAMENBR_TWOSIDED", "all_ds_sample_around_TADs_sameNbr.Rdata")
stopifnot(file.exists(all_ds_sample_around_TADs_sameKb_file))
stopifnot(file.exists(all_ds_sample_around_TADs_fixKb_file))
stopifnot(file.exists(all_ds_sample_around_TADs_sameNbr_file))
all_ds_sample_around_TADs_sameKb <- eval(parse(text = load(all_ds_sample_around_TADs_sameKb_file)))
all_ds_sample_around_TADs_fixKb <- eval(parse(text = load(all_ds_sample_around_TADs_fixKb_file)))
all_ds_sample_around_TADs_sameNbr <- eval(parse(text = load(all_ds_sample_around_TADs_sameNbr_file)))


all_ds <- names(all_ds_sample_around_TADs_sameKb)
stopifnot(setequal(names(all_ds_sample_around_TADs_fixKb), all_ds))
stopifnot(setequal(names(all_ds_sample_around_TADs_sameNbr), all_ds))

ds=all_ds[1]

cat("... prepare all_data_DT_sameKb \n")
all_data_DT_sameKb <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  ds_data <- all_ds_sample_around_TADs_sameKb[[paste0(ds)]]
  ds_data_sub <- lapply(ds_data, function(tad_data) {
    to_keep <- ! names(tad_data) %in% c("tad_genes", "genes", "genes_right", "genes_left")
    tad_data[to_keep]
  })
  ds_data_DT <- do.call(rbind, lapply(ds_data_sub, function(x) data.frame(x)))
  ds_data_DT$region <- rownames(ds_data_DT)
  rownames(ds_data_DT) <- NULL
  ds_data_DT$dataset <- ds
  ds_data_DT
}
all_data_DT_sameKb$sampType <- "sameKb"
colnames(all_data_DT_sameKb)
# [1] "bpAroundTAD"   "nGenes"        "minDist"       "maxDist"       "nGenes_left"  
# [6] "minDist_left"  "maxDist_left"  "nGenes_right"  "minDist_right" "maxDist_right"
# [11] "region"        "dataset"       "sampType"     

cat("... prepare all_data_DT_fixKb \n")
all_data_DT_fixKb <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  ds_data <- all_ds_sample_around_TADs_fixKb[[paste0(ds)]]
  ds_data_sub <- lapply(ds_data, function(tad_data) {
    to_keep <- ! names(tad_data) %in% c("tad_genes", "genes", "genes_right", "genes_left")
    tad_data[to_keep]
  })
  ds_data_DT <- do.call(rbind, lapply(ds_data_sub, function(x) data.frame(x)))
  ds_data_DT$region <- rownames(ds_data_DT)
  rownames(ds_data_DT) <- NULL
  ds_data_DT$dataset <- ds
  ds_data_DT
}
all_data_DT_fixKb$sampType <- "fixKb"
# colnames(all_data_DT_fixKb)
# [1] "nGenes"        "minDist"       "maxDist"       "nGenes_left"   "minDist_left" 
# [6] "maxDist_left"  "nGenes_right"  "minDist_right" "maxDist_right" "region"       
# [11] "dataset"       "sampType"     

cat("... prepare all_data_DT_sameNbr \n")
all_data_DT_sameNbr <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  ds_data <- all_ds_sample_around_TADs_sameNbr[[paste0(ds)]]
  ds_data_sub <- lapply(ds_data, function(tad_data) {
    to_keep <- ! names(tad_data) %in% c("tad_genes", "genes", "genes_right", "genes_left")
    tad_data[to_keep]
  })
  ds_data_DT <- do.call(rbind, lapply(ds_data_sub, function(x) data.frame(x)))
  ds_data_DT$region <- rownames(ds_data_DT)
  rownames(ds_data_DT) <- NULL
  ds_data_DT$dataset <- ds
  ds_data_DT
}
all_data_DT_sameNbr$sampType <- "sameNbr"
# colnames(all_data_DT_sameNbr)
# [1] "nGenes"        "minDist"       "maxDist"       "nGenes_left"   "minDist_left" 
# [6] "maxDist_left"  "nGenes_right"  "minDist_right" "maxDist_right" "region"       
# [11] "dataset"       "sampType"  

commonCols <- Reduce(intersect, list(colnames(all_data_DT_sameNbr),colnames(all_data_DT_fixKb),colnames(all_data_DT_sameKb))) 
# [1] "nGenes"        "minDist"       "maxDist"       "nGenes_left"   "minDist_left" 
# [6] "maxDist_left"  "nGenes_right"  "minDist_right" "maxDist_right" "region"       
# [11] "dataset"       "sampType"     

cat("... prepare all_data_DT \n")
all_data_DT <- rbind(all_data_DT_fixKb[, commonCols] , 
                     rbind(all_data_DT_sameNbr[, commonCols], all_data_DT_sameKb[, commonCols]))



all_vars <- commonCols[! commonCols %in% c("dataset", "region", "sampType")]

plot_var = all_vars[1]

foo <- foreach(plot_var = all_vars) %dopar% {
  cat(paste0("... start plotting for ", plot_var, "\n"))
  stopifnot(plot_var %in% colnames(all_data_DT))
  plotList <- list(
    all_data_DT[all_data_DT$sampType == "fixKb", paste0(plot_var)],
    all_data_DT[all_data_DT$sampType == "sameKb", paste0(plot_var)],
    all_data_DT[all_data_DT$sampType == "sameNbr", paste0(plot_var)]
  )
  names(plotList) <- c("fixKb", "sameKb", "sameNbr")
  
  outFile <- file.path(outFold, paste0(plot_var, "_multiDens.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
  plot_multiDens(
      plotList,
      my_xlab = paste0(plot_var),
     plotTit = paste0(plot_var)
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))

  plotList_log10 <- lapply(plotList, log10)
  outFile <- file.path(outFold, paste0(plot_var, "_multiDens_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
  plot_multiDens(
    plotList_log10,
    my_xlab = paste0(plot_var),
    plotTit = paste0(plot_var)
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
    
}

all_data_DT$dataset_sampType <- paste0(all_data_DT$dataset, "-", all_data_DT$sampType)

cat("... prepare genesCount_byDataSamp_DT \n")
genesCount_byDataSamp_DT <- do.call(rbind, by(all_data_DT, all_data_DT$dataset_sampType, function(subDT) {
  data.frame(
    # totTADs = nrow(subDT),
    no_genes = sum(subDT$nGenes == 0),
    no_genes_left = sum(subDT$nGenes_right == 0),
    no_genes_right = sum(subDT$nGenes_left == 0)
  )
}))
# boxplot(genesCount_byDataSamp_DT)

genesCount_byDataSamp_DT$dataset_sampType <- rownames(genesCount_byDataSamp_DT)
rownames(genesCount_byDataSamp_DT) <- NULL
genesCount_byDataSamp_DT$dataset <- gsub("^(.+)-.+$", "\\1", genesCount_byDataSamp_DT$dataset_sampType)
genesCount_byDataSamp_DT$sampType <- gsub("^.+-(.+)$", "\\1", genesCount_byDataSamp_DT$dataset_sampType)
genesCount_byDataSamp_DT$dataset_sampType <- NULL

outFile <- file.path(outFold, "genesCount_byDataSamp_DT.Rdata")
save(genesCount_byDataSamp_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

plot_genesCount_byDataSamp_DT <- melt(genesCount_byDataSamp_DT, id =c("dataset", "sampType"))

nDS <- length(unique(plot_genesCount_byDataSamp_DT$dataset))

outFile <- file.path(outFold, paste0("nGenes_bySampType_boxplot.", plotType))
p <- ggplot(plot_genesCount_byDataSamp_DT, aes(x = variable, y = value, fill = sampType)) + 
  geom_boxplot() +
  scale_fill_ipsum() +
  scale_color_ipsum() +
  theme_ipsum_rc(grid="XY", axis_title_just="c", axis_title_size=14) + 
  labs(x="", 
       y=paste0(),
       title=paste0("# of genes"),
       # caption="Brought to you by the letter 'g'",
       subtitle=paste0("nDS = ", nDS),
       colour = "",
       fill = ""
  ) 
ggsave(p, filename = outFile, height = myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))





txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))













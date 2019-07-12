options(scipen=100)

script_name <- "look_sampling_corr.R"
cat("... start ", script_name, "\n")
startTime <- Sys.time()


require(reshape2)
require(foreach)
require(ggplot2)
require(hrbrthemes)
require(ggthemes)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

fixKbSize <- 1000000

axisCex <- 1.2

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight

outFold <- "LOOK_SAMPLING_CORR"

all_ds_corr_around_TADs_sameKb_file <- file.path("COEXPR_AROUND_TADS", "sameKb", "all_ds_around_TADs_corr.Rdata")
all_ds_corr_around_TADs_fixKb_file <- file.path("COEXPR_AROUND_TADS", "fixKb", fixKbSize, "all_ds_around_TADs_corr.Rdata")
all_ds_corr_around_TADs_sameNbr_file <- file.path("COEXPR_AROUND_TADS", "sameNbr", "all_ds_around_TADs_corr.Rdata")
stopifnot(file.exists(all_ds_corr_around_TADs_sameKb_file))
stopifnot(file.exists(all_ds_corr_around_TADs_fixKb_file))
stopifnot(file.exists(all_ds_corr_around_TADs_sameNbr_file))
all_ds_corr_around_TADs_sameKb <- eval(parse(text = load(all_ds_corr_around_TADs_sameKb_file)))
all_ds_corr_around_TADs_fixKb <- eval(parse(text = load(all_ds_corr_around_TADs_fixKb_file)))
all_ds_corr_around_TADs_sameNbr <- eval(parse(text = load(all_ds_corr_around_TADs_sameNbr_file)))


all_ds <- names(all_ds_corr_around_TADs_sameKb)
stopifnot(setequal(names(all_ds_corr_around_TADs_fixKb), all_ds))
stopifnot(setequal(names(all_ds_corr_around_TADs_sameNbr), all_ds))

ds=all_ds[1]


all_data_DT_sameKb <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  ds_data <- all_ds_corr_around_TADs_sameKb[[paste0(ds)]]
  ds_data_DT <- do.call(rbind, lapply(ds_data, function(x) data.frame(x)))
  ds_data_DT$region <- rownames(ds_data_DT)
  rownames(ds_data_DT) <- NULL
  ds_data_DT$dataset <- ds
  ds_data_DT
}
all_data_DT_sameKb$sampType <- "sameKb"


all_data_DT_fixKb <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  ds_data <- all_ds_corr_around_TADs_fixKb[[paste0(ds)]]
  ds_data_DT <- do.call(rbind, lapply(ds_data, function(x) data.frame(x)))
  ds_data_DT$region <- rownames(ds_data_DT)
  rownames(ds_data_DT) <- NULL
  ds_data_DT$dataset <- ds
  ds_data_DT
}
all_data_DT_fixKb$sampType <- "fixKb"


all_data_DT_sameNbr <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  ds_data <- all_ds_corr_around_TADs_sameNbr[[paste0(ds)]]
  ds_data_DT <- do.call(rbind, lapply(ds_data, function(x) data.frame(x)))
  ds_data_DT$region <- rownames(ds_data_DT)
  rownames(ds_data_DT) <- NULL
  ds_data_DT$dataset <- ds
  ds_data_DT
}
all_data_DT_sameNbr$sampType <- "sameNbr"


commonCols <- Reduce(intersect, list(colnames(all_data_DT_sameNbr),colnames(all_data_DT_fixKb),colnames(all_data_DT_sameKb))) 


all_data_DT <- rbind(all_data_DT_fixKb[, commonCols] , 
                     rbind(all_data_DT_sameNbr[, commonCols], all_data_DT_sameKb[, commonCols]))



all_vars <- commonCols[! commonCols %in% c("dataset", "region", "sampType")]

plot_var = all_vars[1]

foo <- foreach(plot_var = all_vars) %dopar% {
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
  

}

all_data_DT$dataset_sampType <- paste0(all_data_DT$dataset, "-", all_data_DT$sampType)

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

plot_genesCount_byDataSamp_DT <- melt(genesCount_byDataSamp_DT, id =c("dataset", "sampType"))

nDS <- length(unique(plot_genesCount_byDataSamp_DT$dataset))

outFile <- file.path(outFold, paste0("nGenes_bySampType_boxplot.", plotType))
ggplot(plot_genesCount_byDataSamp_DT, aes(x = variable, y = value, fill = sampType)) + 
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
ggsave(p, filename = outFile, height = myHeight, width = myWidth)
cat(paste0("... written: ", outFile, "\n"))


pipOutFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "OUTPUT_FOLDER")

all_corr_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanCorr_TAD.Rdata", full.names = TRUE)
stopifnot(length(all_corr_files) > 0)


corr_file = all_corr_files[1]
mC_DT <- foreach(corr_file = all_corr_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(corr_file)
  stopifnot(file.exists(curr_file))
  tad_mC <- eval(parse(text = load(curr_file)))
  exprds <- basename(dirname(dirname(corr_file)))
  hicds <- basename(dirname(dirname(dirname(corr_file))))
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mC),
    meanCorr_obs = as.numeric(tad_mC),
    stringsAsFactors = FALSE
  )
}

plot_var = all_vars[1]

plot_var = "meanCorr"
stopifnot(plot_var %in% colnames(all_data_DT))
plotList <- list(
  all_data_DT[all_data_DT$sampType == "fixKb", paste0(plot_var)],
  all_data_DT[all_data_DT$sampType == "sameKb", paste0(plot_var)],
  all_data_DT[all_data_DT$sampType == "sameNbr", paste0(plot_var)],
  mC_DT$meanCorr
)
names(plotList) <- c("fixKb", "sameKb", "sameNbr", "obs.")

outFile <- file.path(outFold, paste0(plot_var, "_withObs_multiDens.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
plot_multiDens(
  plotList,
  my_xlab = paste0(plot_var),
  plotTit = paste0(plot_var)
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

 

all_data_DT$hicds <- dirname(all_data_DT$dataset)
all_data_DT$exprds <- basename(all_data_DT$dataset)


all_data_DT_withObs <- merge(all_data_DT, mC_DT, by = c("hicds", "exprds", "region"), all=TRUE)
stopifnot(nrow(mC_DT) == nrow(all_data_DT))

foo <- foreach(sampType = c("fixKb", "sameKb", "sameNbr")) %dopar% {
  
  
  plotDT <- all_data_DT_withObs[all_data_DT_withObs$sampType == sampType,]
  
  xvar <- "meanCorr_obs"
  myx <- plotDT[,paste0(xvar)]
  
  for(plot_var in all_vars) {
    
    myy <- plotDT[,paste0(plot_var)]
    
    outFile <- file.path(outFolder, paste0(yvar, "_vs_", xvar, "_densplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myHeight))
    densplot(
      x = myx,
      xlab = paste0(xvar),
      y = myy,
      ylab = paste0(plot_var),
      main = paste0(plot_var, " vs. ", xvar),
      cex = 0.7,
      cex.lab = axisCex,
      cex.axis = axisCex
    )
    mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
    addCorr(x=myx, legPos="topleft",
            y=myy, bty='n')
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
  }
  
}





txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))





















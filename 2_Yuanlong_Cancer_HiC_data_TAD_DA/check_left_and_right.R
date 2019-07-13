
# Rscript check_left_and_right.R

options(scipen=100)


script_name <- "check_left_and_right.R"
cat("... start ", script_name, "\n")
startTime <- Sys.time()

source("utils_fct.R")


require(foreach)
require(doMC)

registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
cexPlot <- 1.2

all_sampTypes <- commandArgs(trailingOnly = TRUE)
all_sampTypes = c("sameKb", "fixKb", "sameNbr")

pipFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")

script0_name <- "0_prepGeneData"

buildAll <- FALSE

outFold <- "CHECK_LEFT_AND_RIGHT"
dir.create(outFold, recursive = TRUE)


if(buildAll){
  
  all_dt <- foreach(samp_type = all_sampTypes, .combine='rbind') %do% {
    
    if(samp_type == "fixKb") {
      fixKbSize <- 1000000
    } else{
      fixKbSize <- ""
    }
    check_kb <- ifelse(samp_type == "fixKb", fixKbSize,
                       ifelse(samp_type == "sameNbr", NA, "sameKb"))
    samp_file <- file.path(paste0("CREATE_SAMPLE_AROUND_TADS_", toupper(samp_type), "_TWOSIDED"), 
                           fixKbSize, paste0("all_ds_sample_around_TADs_", samp_type, ".Rdata"))
    stopifnot(file.exists(samp_file))
    all_ds_sample_around <- eval(parse(text = load(samp_file)))
    
    all_ds <- names(all_ds_sample_around)
    
    ds="Panc1_rep12_40kb/TCGApaad_wt_mutKRAS"
    # all_ds=all_ds[1:2]
    
    ds_dt <- foreach(ds = all_ds, .combine='rbind') %do% {
    
      hicds <- dirname(ds)
      exprds <- basename(ds)
      
      stopifnot(dir.exists(file.path(pipFolder, hicds)))
      
      geneList_file <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      geneList <- eval(parse(text = load(geneList_file)))
      
      regionList_file <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
      stopifnot(file.exists(regionList_file))
      regionList <- eval(parse(text = load(regionList_file)))
      
      g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      
      ### RETRIEVE THE TAD POSITIONS
      tadposFile <- file.path(pipFolder, hicds, "genes2tad", "all_assigned_regions.txt")
      stopifnot(file.exists(tadposFile))
      tadpos_DT <- read.delim(tadposFile, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
      stopifnot(is.numeric(tadpos_DT$start))
      stopifnot(is.numeric(tadpos_DT$end))
      tadpos_DT <- tadpos_DT[grepl("_TAD", tadpos_DT$region),,drop=FALSE] 
      
      all_regs <- names(all_ds_sample_around[[paste0(ds)]])
      reg <- all_regs[1]
      # all_regs=all_regs[1:10]
      
      region_dt <- foreach(reg = all_regs, .combine='rbind') %dopar% {
        
  
        cat("... start : ", samp_type, " - ", hicds, " - ", exprds, "  - ", reg, "\n" )
        
        
        stopifnot(reg %in% regionList)
        
        
        
        
        curr_data <- all_ds_sample_around[[paste0(ds)]][[paste0(reg)]]
        
        nGenes <- sum(g2t_DT$region == reg & g2t_DT$entrezID %in% geneList)
        stopifnot(nGenes == length(curr_data$tad_genes))
        
        
        curr_chromo <- gsub("(chr.+)_TAD.+", "\\1", reg)
        
        if(curr_data$nGenes > 0) {
          
          cat("... checkGenes : ", samp_type, " - ", hicds, " - ", exprds, "  - ", reg, "\n" )
          
          checkGenes(genes=curr_data$genes,
                     side="both", region_name=reg, gene2tad_dt=g2t_DT, tadpos_dt=tadpos_DT,
                     checkDistKb = check_kb)
          
          
        }
        if(curr_data$nGenes_left > 0) {
          
          cat("... checkGenes_left : ", samp_type, " - ", hicds, " - ", exprds, "  - ", reg, "\n" )
          
          checkGenes(genes=curr_data$genes_left,
                     side="left", region_name=reg, gene2tad_dt=g2t_DT, tadpos_dt=tadpos_DT,
                     checkDistKb = check_kb)
          
          
        }
        if(curr_data$nGenes_right > 0) {
          
          cat("... checkGenes_right : ", samp_type, " - ", hicds, " - ", exprds, "  - ", reg, "\n" )
          
          checkGenes(genes=curr_data$genes_right,
                     side="right", region_name=reg, gene2tad_dt=g2t_DT, tadpos_dt=tadpos_DT,
                     checkDistKb = check_kb)
          
          
          
        }
        # if(samp_type == "sameKb") {
          tadSize <- (tadpos_DT$end[tadpos_DT$region == reg] - tadpos_DT$start[tadpos_DT$region == reg]) + 1
          stopifnot(length(tadSize) == 1)
          out_dt <- data.frame(
            hicds = hicds,
            exprds = exprds,
            sampType = samp_type,
            tad_size = tadSize,
            max_dist = curr_data$maxDist,
            max_dist_left = curr_data$maxDist_left,
            max_dist_right = curr_data$maxDist_right,
            nGenes_tad = length(curr_data$tad_genes),
            nGenes_sample = curr_data$nGenes,
            nGenes_sample_left = curr_data$nGenes_left,
            nGenes_sample_right = curr_data$nGenes_right,
            stringsAsFactors = FALSE
          )
          return(out_dt)
        #} else {
        #   return(NULL)
        # }
      } # end iterating over regions
      region_dt
        
    } # end-iterating over DS
    ds_dt
  
  }  # end-iterating over kind of sampling
  
  
  outFile <- file.path(outFold, "all_dt.Rdata")
  save(all_dt, file = outFile)
  cat(paste0("... written: ", outFile))

# look at max dist and tad size
} else {
  outFile <- file.path(outFold, "all_dt.Rdata")
  all_dt <- eval(parse(text = load(outFile)))
}


#############################################################################################################################
#############################################################################################################################


for(samp_type in all_sampTypes) {

  sub_all_dt <- all_dt[all_dt$sampType == samp_type,]
  
  
  all_y <- c("max_dist", "max_dist_right", "max_dist_left")
  x_var <- "tad_size"
  y_var = all_y[1]
  
  
  
  foo <- foreach(y_var = all_y) %dopar% {
    
    outFile <- file.path(outFold, paste0(y_var, "_vs_", x_var, "_", samp_type, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = sub_all_dt[,paste0(x_var)],
      y = sub_all_dt[,paste0(y_var)],
      cex.axis = cexPlot,
      cex.lab = cexPlot,
      xlab = paste0(x_var),
      ylab = paste0(y_var),
      main = paste0(y_var, " vs. ", x_var)
    )
    mtext(side=3, text = paste0(samp_type), font=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
   
     
    outFile <- file.path(outFold, paste0(y_var, "_vs_", x_var, "_log10_", samp_type, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = log10(sub_all_dt[,paste0(x_var)]),
      y = log10(sub_all_dt[,paste0(y_var)]),
      cex.axis = cexPlot,
      cex.lab = cexPlot,
      xlab = paste0(x_var, " [log10]"),
      ylab = paste0(y_var, " [log10]"),
      main = paste0(y_var, " vs. ", x_var)
    )
    mtext(side=3, text = paste0(samp_type), font=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  }

  all_y <- c("nGenes_sample", "nGenes_sample_left", "nGenes_sample_right")
  x_var <- "nGenes_tad"
  
  foo <- foreach(y_var = all_y) %dopar% {
    
    outFile <- file.path(outFold, paste0(y_var, "_vs_", x_var,"_", samp_type, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = sub_all_dt[,paste0(x_var)],
      y = sub_all_dt[,paste0(y_var)],
      cex.axis = cexPlot,
      cex.lab = cexPlot,
      xlab = paste0(x_var),
      ylab = paste0(y_var),
      main = paste0(y_var, " vs. ", x_var)
    )
    mtext(side=3, text = paste0(samp_type), font=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    outFile <- file.path(outFold, paste0(y_var, "_vs_", x_var, "_log10_", samp_type, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = log10(sub_all_dt[,paste0(x_var)]),
      y = log10(sub_all_dt[,paste0(y_var)]),
      cex.axis = cexPlot,
      cex.lab = cexPlot,
      xlab = paste0(x_var, " [log10]"),
      ylab = paste0(y_var, " [log10]"),
      main = paste0(y_var, " vs. ", x_var)
    )
    mtext(side=3, text = paste0(samp_type), font=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  
}



}


#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
# all_ds_corr_around_TADs_sameKb_file <- file.path("COEXPR_AROUND_TADS", "sameKb", "all_ds_around_TADs_corr.Rdata")
# all_ds_corr_around_TADs_fixKb_file <- file.path("COEXPR_AROUND_TADS", "fixKb", fixKbSize, "all_ds_around_TADs_corr.Rdata")
# all_ds_corr_around_TADs_sameNbr_file <- file.path("COEXPR_AROUND_TADS", "sameNbr", "all_ds_around_TADs_corr.Rdata")
# stopifnot(file.exists(all_ds_corr_around_TADs_sameKb_file))
# stopifnot(file.exists(all_ds_corr_around_TADs_fixKb_file))
# stopifnot(file.exists(all_ds_corr_around_TADs_sameNbr_file))
# all_ds_corr_around_TADs_sameKb <- eval(parse(text = load(all_ds_corr_around_TADs_sameKb_file)))
# all_ds_corr_around_TADs_fixKb <- eval(parse(text = load(all_ds_corr_around_TADs_fixKb_file)))
# all_ds_corr_around_TADs_sameNbr <- eval(parse(text = load(all_ds_corr_around_TADs_sameNbr_file)))
# 
# ### sameKb - meanCorr
# all_meanCorr_sameKb <- unlist(lapply(all_ds_corr_around_TADs_sameKb, function(sub_data){
#   lapply(sub_data, function(x) x[["meanCorr"]])
# }))
# all_meanCorr_sameKb <- na.omit(all_meanCorr_sameKb)
# 
# ### sameKb - meanCorrRight
# all_meanCorrRight_sameKb <- unlist(lapply(all_ds_corr_around_TADs_sameKb, function(sub_data){
#   lapply(sub_data, function(x) x[["meanCorr_right"]])
# }))
# all_meanCorrRight_sameKb <- na.omit(all_meanCorrRight_sameKb)
# 
# ### sameKb - meanCorr
# all_meanCorrLeft_sameKb <- unlist(lapply(all_ds_corr_around_TADs_sameKb, function(sub_data){
#   lapply(sub_data, function(x) x[["meanCorr_left"]])
# }))
# all_meanCorrLeft_sameKb <- na.omit(all_meanCorrLeft_sameKb)
# 
# 
# all_empPvalsLeft_sameKb <- sapply(all_obs_corr, function(x) {
#   (sum(all_meanCorrLeft_sameKb >= x) + 1)/length(all_meanCorrLeft_sameKb)
# })
# 

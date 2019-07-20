# Rscript create_sample_alternative.R 10
# Rscript create_sample_alternative.R 1000

# for each TAD => sample "n" random TADs on the same chromo (same number of adjacent genes)
# => same as going on a chromo, for each TAD size sample nTAD_of_this_size * n

set.seed(20072019)

script_name <- "create_sample_alternative.R"

cat("... start ", script_name, "\n")

startTime <- Sys.time()

options(scipen=100)

plotType <- "png"
myHeight <- 400
myWidth <- 600

SSHFS <- FALSE
nCpu <- ifelse(SSHFS, 2, 40)


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(nCpu)

nSampleByTAD=10
nSampleByTAD <- commandArgs(trailingOnly = TRUE)[1]
nSampleByTAD <- as.numeric(nSampleByTAD)
stopifnot(!is.na(nSampleByTAD))


outFolder <- file.path("CREATE_SAMPLE_ALTERNATIVE", nSampleByTAD)
dir.create(outFolder, recursive = TRUE)

pipfolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")

pipOutFolder <- file.path(pipfolder, "PIPELINE", "OUTPUT_FOLDER")

stopifnot(dir.exists(pipOutFolder))

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(length(all_hicexpr_ds) > 0)
stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))

ds="ENCSR862OGI_RPMI-7951_40kb/TCGAskcm_lowInf_highInf"
ds=all_hicexpr_ds[1]

ds="ENCSR401TBQ_Caki2_40kb/TCGAkich_norm_kich"
reg = "chr9_TAD80"

#all_hicexpr_ds = all_hicexpr_ds[1:3]
# all_hicexpr_ds=all_hicexpr_ds[1]

all_ds_sample_byTAD <-   foreach(ds = all_hicexpr_ds) %do% {
  hicds <- file.path(dirname(ds))
  exprds <- basename(ds)
  stopifnot(dir.exists(file.path(pipfolder, hicds)))
  
  cat("... start dataset: ", exprds, "\n")
  
  dsPipOutDir <- file.path(pipOutFolder, ds)
  stopifnot(dir.exists(dsPipOutDir))
  
  ### RETRIEVE THE GENE2TAD ASSIGNMENT
  g2tFile <- file.path(pipfolder, hicds, "genes2tad", "all_genes_positions.txt")
  stopifnot(file.exists(g2tFile))
  g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
  g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
  
  ### RETRIEVE THE TAD POSITIONS
  tadposFile <- file.path(pipfolder, hicds, "genes2tad", "all_assigned_regions.txt")
  stopifnot(file.exists(tadposFile))
  tadpos_DT <- read.delim(tadposFile, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
  stopifnot(is.numeric(tadpos_DT$start))
  stopifnot(is.numeric(tadpos_DT$end))
  tadpos_DT <- tadpos_DT[grepl("_TAD", tadpos_DT$region),,drop=FALSE] 
  
  
  ### KEEP ONLY THE TADs USED IN THE PIPELINE
  script0_name <- "0_prepGeneData"
  stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
  tadListFile <- file.path(dsPipOutDir, script0_name, "pipeline_regionList.Rdata")
  stopifnot(file.exists(tadListFile))
  pipeline_tadList <- eval(parse(text = load(tadListFile))) # not adjusted
  stopifnot(pipeline_tadList %in% tadpos_DT$region)
  tadpos_DT <- tadpos_DT[tadpos_DT$region %in% pipeline_tadList,]
  
  stopifnot(!duplicated(pipeline_tadList))
  
  ### RETRIEVE THE GENES USED IN THE PIPELINE - script0
  script0_name <- "0_prepGeneData"
  stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
  geneListFile <- file.path(dsPipOutDir, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(geneListFile))
  pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted
  stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
  
  stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
  # stopifnot(names(pipeline_geneList) %in% g2t_DT$entrezID) -> FALSE
  g2t_DT <- g2t_DT[as.character(g2t_DT$entrezID) %in% as.character(pipeline_geneList),,drop=FALSE]
  stopifnot(length(pipeline_geneList) == nrow(g2t_DT))
  stopifnot(g2t_DT$entrezID %in% pipeline_geneList)
  g2t_DT$chromo <- as.character(g2t_DT$chromo)
  
  stopifnot(g2t_DT$entrezID %in% pipeline_geneList)
  stopifnot(grepl("TAD", g2t_DT$region))
  
  tadpos_DT$mid_pos <- (tadpos_DT$start+tadpos_DT$end)/2
  g2t_DT$mid_pos <- (g2t_DT$start+g2t_DT$end)/2
  
  ### !!! FILTERED g2t BY GENELIST !!!!!!!!!!!!!!
  
  all_chromo <- unique(g2t_DT$chromo)
  
  chromo=all_chromo[1]
  all_chr_sampling <- foreach(chromo = all_chromo) %do% {
    
    cat("... start dataset: ", exprds, " - ", chromo, "\n")
    
    
    chr_g2t_DT <- g2t_DT[g2t_DT$chromo == chromo,]
    stopifnot(nrow(chr_g2t_DT) > 0)
    
    # get how many by tad size
    genesByTAD <- setNames(as.numeric(table(chr_g2t_DT$region)), as.character(names(table(chr_g2t_DT$region))))
    
    TADsizeDist <- setNames(as.numeric(table(genesByTAD)), as.character(names(table(genesByTAD))))
    stopifnot(length(TADsizeDist) > 0)
    
    # do the sampling for each size
    i=1
    all_size_sampling <- foreach(i = seq_along(TADsizeDist)) %dopar% {
      
      cat("... start dataset: ", exprds, " - ", chromo, " - size=", names(TADsizeDist[i]), "\n")
      
      nGenesToSample <- as.numeric(as.character(names(TADsizeDist)[i]))
      
      nSampling <- nSampleByTAD * as.numeric(TADsizeDist[i])
      stopifnot(!is.na(nSampling))
      stopifnot(nSampling > 0)
      
      maxIdx <- nrow(chr_g2t_DT)-nGenesToSample+1
      if(maxIdx < nSampling) {
        maxSampling <- maxIdx
      } else {
        maxSampling <- nSampling
      }
      
      all_idxs <- sample(c(1:maxIdx), maxSampling, replace=FALSE) # if a DT of 40 rows and 3 genes to sample, go up to 38 
      stopifnot(!duplicated(all_idxs))
      
      k_idx = all_idxs[1]
      size_sampling <- foreach(k_idx = all_idxs) %do% {
        
        
        
        stopifnot(is.numeric(k_idx))
        
        last_idx <- k_idx+nGenesToSample-1
        stopifnot(last_idx <= nrow(chr_g2t_DT))
        
        
        selected_DT <- chr_g2t_DT[(k_idx:last_idx),]
        
        if(length(unique(selected_DT$region)) == 1) return(NULL)
        
        selected_genes <- selected_DT$entrezID
        
      
        stopifnot(selected_genes %in% pipeline_geneList)
        selected_genes
        
        
      } # end-foreach iterating over the sampling to perform for the current TAD size of the current chromo
      names(size_sampling) <- paste0("samp_", 1:length(all_idxs))
      size_sampling
      
    } # end-foreach iterating over the size of the TADs for the current chromo
    names(all_size_sampling) <- paste0(TADsizeDist, "_TADsOfSize_", names(TADsizeDist))
    all_size_sampling
    
    
    
  } # end-foreach iterating over the chromo
  names(all_chr_sampling) <- all_chromo
  all_chr_sampling
} # end-foreach iterating over datasets

names(all_ds_sample_byTAD) <- all_hicexpr_ds

  
  

outFile <- file.path(outFolder, "all_ds_sample_byTAD.Rdata")
save(all_ds_sample_byTAD, file=outFile)
cat(paste0("... written: ", outFile, "\n"))



########################################################
########################################################
########################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))









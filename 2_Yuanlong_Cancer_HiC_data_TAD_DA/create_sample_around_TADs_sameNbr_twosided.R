# Rscript create_sample_around_TADs_sameNbr_twosided.R

# for a given number of genes
# go at each BD location
# select "n" genes left and right of the BD
# to get an empirical distribution of the correlation
# between expression of genes separated by the boundary

script_name <- "create_sample_around_TADs_sameNbr_twosided.R"

cat("... start ", script_name, "\n")

startTime <- Sys.time()


plotType <- "png"
myHeight <- 400
myWidth <- 600

SSHFS <- FALSE
nCpu <- ifelse(SSHFS, 2, 40)

buildSampleAroundTADs <- TRUE

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(nCpu)

outFolder <- file.path("CREATE_SAMPLE_AROUND_TADS_SAMENBR_TWOSIDED")
dir.create(outFolder)

pipfolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")

pipOutFolder <- file.path(pipfolder, "PIPELINE", "OUTPUT_FOLDER")

stopifnot(dir.exists(pipOutFolder))

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(length(all_hicexpr_ds) > 0)
stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))

ds="ENCSR862OGI_RPMI-7951_40kb/TCGAskcm_lowInf_highInf"
ds=all_hicexpr_ds[1]
# all_hicexpr_ds=all_hicexpr_ds[1]
#all_hicexpr_ds = all_hicexpr_ds[1:3]

if(buildSampleAroundTADs) {
  
  
  # all_hicexpr_ds=all_hicexpr_ds[1]
  all_ds_sample_around_TADs_sameNbr <-   foreach(ds = all_hicexpr_ds) %do% {
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
        
    # nPos=10
    reg=pipeline_tadList[1]
    all_sample_around_TADs <- foreach(reg = pipeline_tadList) %dopar% {

            
      cat("...... start TAD : \t", reg, "\n")
      
      
      curr_chromo <- as.character(tadpos_DT$chromo[tadpos_DT$region == reg])
      
      curr_start <- (tadpos_DT$start[tadpos_DT$region == reg])
      stopifnot(is.numeric(curr_start))
      
      
      curr_end <- (tadpos_DT$end[tadpos_DT$region == reg])
      stopifnot(is.numeric(curr_end))

      stopifnot(length(curr_chromo) == 1)
      stopifnot(length(curr_start) == 1)
      stopifnot(length(curr_end) == 1)
      
      curr_midPos <- (curr_start+curr_end)/2
      stopifnot(curr_midPos == tadpos_DT$mid_pos[tadpos_DT$region == reg])
      
      reg_genes <- g2t_DT$entrezID[g2t_DT$region == reg]
      stopifnot(length(reg_genes) > 0)
      
      curr_nGenes <- length(reg_genes)
      
      # !!! EXTRACT GENES BASED ON START POSITION RELATIVE TO BD
      # !!! SMALLER THAN / GREATER *OR EQUAL* THAN BD POSITION (smaller not equal otherwise genes could come twice)
      
      curr_g2t <- g2t_DT[g2t_DT$chromo == curr_chromo,,drop=FALSE]
      stopifnot(nrow(curr_g2t) > 0)
      
      stopifnot(is.numeric(curr_g2t$start), is.numeric(curr_g2t$end))
      curr_g2t <- curr_g2t[order(curr_g2t$start, curr_g2t$end),,drop=FALSE]
      
      
      # distance to TAD center
      curr_genesOutsideDT <- curr_g2t[ !curr_g2t$entrezID %in% reg_genes,,drop=FALSE]
      stopifnot(nrow(curr_genesOutsideDT) > 0)
        
      #>>> take the same number of genes, either on left or right
      curr_genesOutsideDT$distToTAD <- abs(curr_genesOutsideDT$mid_pos - curr_midPos)
      curr_genesOutsideDT <- curr_genesOutsideDT[order(curr_genesOutsideDT$distToTAD, decreasing=F),,drop=FALSE]
      
      sample_around_genes <- curr_genesOutsideDT$entrezID[1:curr_nGenes]
      
      all_dist <- curr_genesOutsideDT$distToTAD[1:curr_nGenes]
      stopifnot(!is.na(all_dist))
      
      stopifnot(!is.na(curr_genesOutsideDT))
      
      stopifnot(length(all_dist) == length(sample_around_genes))
      
      #>>> take the same number of genes on the left
      curr_genesOutsideDT_left <- curr_genesOutsideDT[curr_genesOutsideDT$mid_pos < curr_midPos,]  # SMALLER
      curr_nGenes_left <- min(curr_nGenes, nrow(curr_genesOutsideDT_left)) 
      
      if(curr_nGenes_left > 0) {
        sample_around_genes_left <- curr_genesOutsideDT_left$entrezID[1:curr_nGenes_left]
        all_dist_left <- curr_genesOutsideDT_left$distToTAD[1:curr_nGenes_left]  
        stopifnot(!is.na(all_dist_left))
        ### ONLY CONSIDER IN COEXPR GENES USED IN PIPELINE ???
        stopifnot(sample_around_genes_left %in% pipeline_geneList)
        stopifnot(!sample_around_genes_left %in% reg_genes)
        
      } else {
        sample_around_genes_left <- character(0)
        all_dist_left <- c()
      }
      
      
      
      curr_genesOutsideDT_right <- curr_genesOutsideDT[curr_genesOutsideDT$mid_pos >= curr_midPos,]  # BIGGER OR EQUAL
      curr_nGenes_right <- min(curr_nGenes, nrow(curr_genesOutsideDT_right)) 
      
      if(curr_nGenes_right > 0) {
        sample_around_genes_right <- curr_genesOutsideDT_right$entrezID[1:curr_nGenes_right]
        all_dist_right <- curr_genesOutsideDT_right$distToTAD[1:curr_nGenes_right]
        stopifnot(!is.na(all_dist_right))
        
        stopifnot(sample_around_genes_right %in% pipeline_geneList)
        stopifnot(!sample_around_genes_right %in% reg_genes)
        
      } else {
        sample_around_genes_right <- character(0)
        all_dist_right <- c()
        
      }
      
      
      
      stopifnot(length(sample_around_genes) == curr_nGenes)
      stopifnot(length(sample_around_genes_left) <= curr_nGenes)
      stopifnot(length(sample_around_genes_right) <= curr_nGenes)
      
      
      
      stopifnot(length(all_dist) == length(sample_around_genes))
      stopifnot(length(all_dist_right) == length(sample_around_genes_right))
      stopifnot(length(all_dist_left) == length(sample_around_genes_left))
      
      
      ### ONLY CONSIDER IN COEXPR GENES USED IN PIPELINE ???
      stopifnot(sample_around_genes %in% pipeline_geneList)
      
      list(
          tad_genes = reg_genes,
          
           genes = sample_around_genes,
           nGenes = length(sample_around_genes),
           minDist = min(all_dist),
           maxDist = max(all_dist),
          
          
          genes_left = sample_around_genes_left,
          nGenes_left = curr_nGenes_left,
          minDist_left = min(all_dist_left),
          maxDist_left = max(all_dist_left),
          
          genes_right = sample_around_genes_right,
          nGenes_right =  curr_nGenes_right,
          minDist_right = min(all_dist_right),
          maxDist_right = max(all_dist_right)
          
          
          
          )
      
      
      
      } # end foreach-iterating over TADs
    
    names(all_sample_around_TADs) <- pipeline_tadList

  sample_around_TADs_sameNbr <- all_sample_around_TADs
  outFile <- file.path(outFolder, hicds, exprds, "sample_around_TADs_sameNbr.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(sample_around_TADs_sameNbr, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))


    all_sample_around_TADs
  } # end for iterating over ds
  names(all_ds_sample_around_TADs_sameNbr) <- all_hicexpr_ds
  outFile <- file.path(outFolder, "all_ds_sample_around_TADs_sameNbr.Rdata")
  save(all_ds_sample_around_TADs_sameNbr, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  # end-if buildCorrAroundBD
} else {
  outFile <- file.path(outFolder, "all_ds_sample_around_TADs_sameNbr.Rdata")
  all_ds_sample_around_TADs_sameNbr <- eval(parse(text = load(outFile)))
}





########################################################
########################################################
########################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))










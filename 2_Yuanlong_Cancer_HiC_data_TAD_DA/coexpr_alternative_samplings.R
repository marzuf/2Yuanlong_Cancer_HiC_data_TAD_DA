
options(scipen=100)

# Rscript coexpr_alternative_samplings.R 10

script_name <- "coexpr_alternative_samplings.R"

startTime <- Sys.time()

cat("> START ",script_name, "\n")

SSHFS <- FALSE

buildData <- TRUE

setDir <- ifelse(SSHFS, "/media/electron", "") # needed to load the setting file...

require(foreach)
require(reshape2)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 90))
source("utils_fct.R")


nSampled =10
nSampled <- commandArgs(trailingOnly = TRUE)
nSampled <- as.numeric(nSampled)
stopifnot(!is.na(nSampled))


outFolder <- file.path("COEXPR_ALTERNATIVE_SAMPLINGS")
dir.create(outFolder, recursive=TRUE)

inFold <- file.path("CREATE_SAMPLE_ALTERNATIVE", nSampled)
stopifnot(dir.exists(inFold))
inFile <- file.path(inFold, "all_ds_sample_byTAD.Rdata")
stopifnot(file.exists(inFile))
all_ds_sample_data <- eval(parse(text = load(inFile)))

corrMeth <- "pearson"
script0_name <- "0_prepGeneData"

pipFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))




all_ds <- names(all_ds_sample_data)

nDS <- length(all_ds)

cat("... found nDS = ", nDS, "\n")

ds = all_ds[2]


all_ds=all_ds[1]

all_ds_sampling_corr <- foreach(ds = all_ds ) %do% {
  
  cat(paste0("... start DS = ", ds, "\n"))
  
  hicds <- dirname(ds)
  exprds <- basename(ds)
  stopifnot(dir.exists(file.path(pipFolder, hicds)))
  
  geneList_file <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(geneList_file))
  geneList <- eval(parse(text = load(geneList_file)))
  
  
  regionList_file <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
  stopifnot(file.exists(regionList_file))
  regionList <- eval(parse(text = load(regionList_file)))
  
  
  settingFile <- file.path(pipFolder, "PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R"))
  stopifnot(file.exists(settingFile))
  source(settingFile)
  sample1_file <- file.path(setDir, sample1_file)
  sample2_file <- file.path(setDir, sample2_file)
  stopifnot(file.exists(sample1_file))
  stopifnot(file.exists(sample2_file))
  
  cond1_ID <- eval(parse(text = load(sample1_file)))
  cond2_ID <- eval(parse(text = load(sample2_file)))
  
  qqnormDTfile <- file.path(pipOutFolder, hicds, exprds,script0_name, "rna_qqnorm_rnaseqDT.Rdata")
  stopifnot(file.exists(qqnormDTfile))
  qqnormDT <- eval(parse(text = load(qqnormDTfile)))
  
  stopifnot(names(geneList) %in% rownames(qqnormDT))
  stopifnot(setequal(colnames(qqnormDT), c(cond1_ID, cond2_ID)))
  
  norm_rnaseqDT <- qqnormDT[names(geneList),]    # !!! ENSURE THAT THE QQNORM IN THE SAME ORDER AS THE GENELIST !!!
  
  stopifnot(rownames(norm_rnaseqDT) == names(geneList))
  stopifnot(!duplicated(names(geneList)))
  
  ds_sample_data <- all_ds_sample_data[[paste0(ds)]]
  
  all_chr <- names(ds_sample_data)
  
  
  chr = all_chr[7]
  ds_all_corr_data <- foreach(chr = all_chr) %dopar% {
    
    cat(paste0("... start DS = ", ds, " - ", chr, "\n"))
    
    
    chr_sample_data <- ds_sample_data[[paste0(chr)]]
      
    stopifnot(length(chr_sample_data) > 0)
    i=1
    chr_corr_data <- lapply(c(1:length(chr_sample_data)), function(i) {
      
      sample_data <- chr_sample_data[[i]]    
      
      TADsize <- gsub(".+TADsOfSize_(.+)", "\\1", names(chr_sample_data)[i])
      TADsize <- as.numeric(TADsize)
      stopifnot(!is.na(TADsize))
        
      sample_data <- Filter(function(x) !is.null(x), sample_data)
      
      if(length(sample_data) == 0) return(NA)
      sample_genes = sample_data[[1]]
      sample_corr_data <- lapply(sample_data, function(sample_genes) {
        
        stopifnot(length(sample_genes) == TADsize)
        

        
        
        corr_genes <- names(geneList)[geneList %in% sample_genes]      # the names used in the nrom_rnaseqDT
        

        
        
        sub_normDT <- norm_rnaseqDT[c(corr_genes),]
        
        stopifnot(nrow(sub_normDT) == length(corr_genes))
        
        stopifnot(cond1_ID %in% colnames(sub_normDT))
        stopifnot(cond2_ID %in% colnames(sub_normDT))
        
        sub_normDT_cond1 <- sub_normDT[,cond1_ID]
        sub_normDT_cond2 <- sub_normDT[,cond2_ID]
        
        stopifnot(nrow(sub_normDT) == nrow(sub_normDT_cond1))
        stopifnot(nrow(sub_normDT) == nrow(sub_normDT_cond2))
        stopifnot( ncol(sub_normDT_cond1) + ncol(sub_normDT_cond2) == ncol(sub_normDT))
        stopifnot( ncol(sub_normDT_cond1) == length(cond1_ID))
        stopifnot(ncol(sub_normDT_cond2) == length(cond2_ID))
        
        
        
        meanCorr_all <- get_meanCorr_value_alternative(
          exprMatrix = sub_normDT, 
          sample_genes = corr_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
          cormet = corrMeth
        )
        
        meanCorr_cond1 <- get_meanCorr_value_alternative(
          exprMatrix = sub_normDT_cond1, 
          sample_genes = corr_genes, 
          cormet = corrMeth
        )
        
        
        meanCorr_cond2 <- get_meanCorr_value_alternative(
          exprMatrix = sub_normDT_cond2, 
          sample_genes = corr_genes, 
          cormet = corrMeth
        )
        
        
        
        list(meanCorr_all=meanCorr_all,
             meanCorr_cond1 = meanCorr_cond1,
             meanCorr_cond2 = meanCorr_cond2)
        
        
        
        
      }) # end-lapply over sample_data (=list of genes for the current sampling)
      
      outFile <- file.path(outFolder, "sample_corr_data.Rdata")
      save(sample_corr_data, file = outFile)
      cat(paste0("... written: ", outFile, "\n"))
      
      
      return(sample_corr_data)
      
      
    }) # end-lapply over the different TAD size of the current chromo  
    outFile <- file.path(outFolder, "chr_corr_data.Rdata")
    save(chr_corr_data, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    names(chr_corr_data) <- names(chr_sample_data)
    return(chr_corr_data)
  } # end-foreach iterating over the chromo
  names(ds_all_corr_data) <- all_chr
  outFile <- file.path(outFolder, "ds_all_corr_data.Rdata")
  save(ds_all_corr_data, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  ds_all_corr_data
} #end-foreach iterating over the datasets
      
names(all_ds_sampling_corr) <- all_ds
  

outFile <- file.path(outFolder, "all_ds_sampling_corr.Rdata")
save(all_ds_sampling_corr, file = outFile)
cat(paste0("... written: ", outFile, "\n"))



####################################################################################
####################################################################################3
####################################################################################3
txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))


#     
#     if(tad_data$nGenes > 0) {
#       
#       ########## => TAKING THE TADs AND SAMPLING ON BOTH SIDES OF THE TADs
#       sample_genes <- tad_data$genes
#       tad_genes <- tad_data$tad_genes
#       
#       stopifnot(! sample_genes %in% tad_genes)
#       stopifnot(sample_genes %in% geneList)
#       
#       inTAD_genes <- names(geneList)[geneList %in% tad_genes]      # the names used in the nrom_rnaseqDT
#       outTAD_genes <- names(geneList)[geneList %in% sample_genes]  # the names used in the nrom_rnaseqDT
#       
#       nTotGenes <- length(inTAD_genes) + length(outTAD_genes)
#       
#       stopifnot(! inTAD_genes %in% outTAD_genes)
#       stopifnot(! outTAD_genes %in% inTAD_genes)
#       stopifnot(inTAD_genes %in% rownames(norm_rnaseqDT))
#       stopifnot(outTAD_genes %in% rownames(norm_rnaseqDT))
#       
#       sub_normDT <- norm_rnaseqDT[c(inTAD_genes, outTAD_genes),]
#       
#       stopifnot(nrow(sub_normDT) == nTotGenes)
#       stopifnot(rownames(sub_normDT) == c(inTAD_genes, outTAD_genes))
#       stopifnot(cond1_ID %in% colnames(sub_normDT))
#       stopifnot(cond2_ID %in% colnames(sub_normDT))
#       
#       sub_normDT_cond1 <- sub_normDT[,cond1_ID]
#       sub_normDT_cond2 <- sub_normDT[,cond2_ID]
#       
#       stopifnot(nrow(sub_normDT) == nrow(sub_normDT_cond1))
#       stopifnot(nrow(sub_normDT) == nrow(sub_normDT_cond2))
#       stopifnot( ncol(sub_normDT_cond1) + ncol(sub_normDT_cond2) == ncol(sub_normDT))
#       stopifnot( ncol(sub_normDT_cond1) == length(cond1_ID))
#       stopifnot(ncol(sub_normDT_cond2) == length(cond2_ID))
#       
# 
#       
#       meanCorr_all <- get_meanCorr_value(
#                      exprMatrix = sub_normDT, 
#                      inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
#                      outside_genes = outTAD_genes, 
#                      cormet = corrMeth
#                      )
#         
#       meanCorr_cond1 <- get_meanCorr_value(
#         exprMatrix = sub_normDT_cond1, 
#         inside_genes = inTAD_genes, 
#         outside_genes = outTAD_genes, 
#         cormet = corrMeth
#       )
#       
#       
#       meanCorr_cond2 <- get_meanCorr_value(
#         exprMatrix = sub_normDT_cond2, 
#         inside_genes = inTAD_genes, 
#         outside_genes = outTAD_genes, 
#         cormet = corrMeth
#       )
#       
#       
# 
#       ########## => TAKING THE TADs AND SAMPLING ON THE RIGHT ONLY
#       if(tad_data$nGenes_right > 0) {
#         
#         sample_genes_right <- tad_data$genes_right
#         outTAD_genes_right <- names(geneList)[geneList %in% sample_genes_right]  # the names used in the nrom_rnaseqDT
#         
#         stopifnot(! sample_genes_right %in% tad_genes)
#         
#         
#         nTotGenes_right <- length(inTAD_genes) + length(outTAD_genes_right)
#         
#         stopifnot(! inTAD_genes %in% outTAD_genes_right)
#         stopifnot(! outTAD_genes_right %in% inTAD_genes)
# 
#         stopifnot(outTAD_genes_right %in% rownames(norm_rnaseqDT))
#         
#         sub_normDT_right <- norm_rnaseqDT[c(inTAD_genes, outTAD_genes_right),]
#         
#         stopifnot(nrow(sub_normDT_right) == nTotGenes_right)
#         stopifnot(rownames(sub_normDT_right) == c(inTAD_genes, outTAD_genes_right))
#         stopifnot(cond1_ID %in% colnames(sub_normDT_right))
#         stopifnot(cond2_ID %in% colnames(sub_normDT_right))
#         
#         sub_normDT_cond1_right <- sub_normDT_right[,cond1_ID]
#         sub_normDT_cond2_right <- sub_normDT_right[,cond2_ID]
#         
#         stopifnot(nrow(sub_normDT_right) == nrow(sub_normDT_cond1_right))
#         stopifnot(nrow(sub_normDT_right) == nrow(sub_normDT_cond2_right))
#         stopifnot( ncol(sub_normDT_cond1_right) + ncol(sub_normDT_cond2_right) == ncol(sub_normDT_right))
#         stopifnot( ncol(sub_normDT_cond1_right) == length(cond1_ID))
#         stopifnot(ncol(sub_normDT_cond2_right) == length(cond2_ID))
#         
#         
#         meanCorrRight_all <- get_meanCorr_value(
#           exprMatrix = sub_normDT_right, 
#           inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
#           outside_genes = outTAD_genes_right, 
#           cormet = corrMeth
#         )
#         
#         meanCorrRight_cond1 <- get_meanCorr_value(
#           exprMatrix = sub_normDT_cond1_right, 
#           inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
#           outside_genes = outTAD_genes_right, 
#           cormet = corrMeth
#         )
#         
#         meanCorrRight_cond2 <- get_meanCorr_value(
#           exprMatrix = sub_normDT_cond2_right, 
#           inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
#           outside_genes = outTAD_genes_right, 
#           cormet = corrMeth
#         )
#         
#         
#       } else {
#         
#         
#         meanCorrRight_all <- NA
#         meanCorrRight_cond1 <- NA
#         meanCorrRight_cond2 <- NA
#         
#       }
#       if(tad_data$nGenes_left > 0) {
#         
#         sample_genes_left <- tad_data$genes_left
#         
#         stopifnot(! sample_genes_left %in% tad_genes)
#         
#         outTAD_genes_left <- names(geneList)[geneList %in% sample_genes_left]  # the names used in the nrom_rnaseqDT
#         
#         
#         
#         stopifnot(! sample_genes_left %in% tad_genes)
#         
#         
#         nTotGenes_left <- length(inTAD_genes) + length(outTAD_genes_left)
#         
#         stopifnot(! inTAD_genes %in% outTAD_genes_left)
#         stopifnot(! outTAD_genes_left %in% inTAD_genes)
#         
#         stopifnot(outTAD_genes_left %in% rownames(norm_rnaseqDT))
#         
#         sub_normDT_left <- norm_rnaseqDT[c(inTAD_genes, outTAD_genes_left),]
#         
#         stopifnot(nrow(sub_normDT_left) == nTotGenes_left)
#         stopifnot(rownames(sub_normDT_left) == c(inTAD_genes, outTAD_genes_left))
#         stopifnot(cond1_ID %in% colnames(sub_normDT_left))
#         stopifnot(cond2_ID %in% colnames(sub_normDT_left))
#         
#         sub_normDT_cond1_left <- sub_normDT_left[,cond1_ID]
#         sub_normDT_cond2_left <- sub_normDT_left[,cond2_ID]
#         
#         stopifnot(nrow(sub_normDT_left) == nrow(sub_normDT_cond1_left))
#         stopifnot(nrow(sub_normDT_left) == nrow(sub_normDT_cond2_left))
#         stopifnot( ncol(sub_normDT_cond1_left) + ncol(sub_normDT_cond2_left) == ncol(sub_normDT_left))
#         stopifnot( ncol(sub_normDT_cond1_left) == length(cond1_ID))
#         stopifnot(ncol(sub_normDT_cond2_left) == length(cond2_ID))
#         
#         
#         meanCorrLeft_all <- get_meanCorr_value(
#           exprMatrix = sub_normDT_left, 
#           inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
#           outside_genes = outTAD_genes_left, 
#           cormet = corrMeth
#         )
#         
#         meanCorrLeft_cond1 <- get_meanCorr_value(
#           exprMatrix = sub_normDT_cond1_left, 
#           inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
#           outside_genes = outTAD_genes_left, 
#           cormet = corrMeth
#         )
#         
#         meanCorrLeft_cond2 <- get_meanCorr_value(
#           exprMatrix = sub_normDT_cond2_left, 
#           inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
#           outside_genes = outTAD_genes_left, 
#           cormet = corrMeth
#         )
#         
#         
#       } else {
#         
#         
#         meanCorrLeft_all <- NA
#         meanCorrLeft_cond1 <- NA
#         meanCorrLeft_cond2 <- NA
#         
#       }
#       
#       
#       
#     } else {
#       meanCorr_all <- NA
#       meanCorr_cond1 <- NA
#       meanCorr_cond2 <- NA
#       
#       meanCorrRight_all <- NA
#       meanCorrRight_cond1 <- NA
#       meanCorrRight_cond2 <- NA
#       
#       
#       meanCorrLeft_all <- NA
#       meanCorrLeft_cond1 <- NA
#       meanCorrLeft_cond2 <- NA
#       
#     }
#     
#     list(
#       nGenes = tad_data$nGenes,
#       meanCorr = meanCorr_all,
#       meanCorr_cond1 = meanCorr_cond1,
#       meanCorr_cond2 = meanCorr_cond2,
#       nGenes_right = tad_data$nGenes_right,
#       meanCorr_right = meanCorrRight_all,
#       meanCorr_cond1_right = meanCorrRight_cond1,
#       meanCorr_cond2_right = meanCorrRight_cond2,
#       nGenes_left = tad_data$nGenes_left,
#       meanCorr_left = meanCorrLeft_all,
#       meanCorr_cond1_left = meanCorrLeft_cond1,
#       meanCorr_cond2_left = meanCorrLeft_cond2
#     )
#     
#   } # end iterating over all TADs for the current dataset
#   names(ds_all_corr_data) <- all_regs
#   ds_all_corr_data
# } # end iterating over all DS
# names(all_ds_around_TADs_corr) <- all_ds 
#   
# 


# Rscript coexpr_around_TADs.R fixKb 1000000
# Rscript coexpr_around_TADs.R sameKb
# Rscript coexpr_around_TADs.R sameNbr

script_name <- "coexpr_around_TADs.R"

startTime <- Sys.time()

cat("> START coexpr_around_TADs.R \n")

SSHFS <- FALSE

buildData <- TRUE

setDir <- ifelse(SSHFS, "/media/electron", "") # needed to load the setting file...

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 90))


withDiag <- FALSE

# CREATE_SAMPLE_AROUND_TADS_FIXKB  CREATE_SAMPLE_AROUND_TADS_SAMEKB   CREATE_SAMPLE_AROUND_TADS_SAMENBR
# args <- c("sameKb")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(args[1] %in% c("fixKb", "sameKb", "sameNbr"))
if(args[1] == "fixKb") {
  stopifnot(length(args) == 2)
  stopifnot(!is.na(as.numeric(as.character(args[2]))))
  windowSizeBp <- args[2]
} else {
  args[2] <- ""
}

options(scipen=100)

outFolder <- file.path("COEXPR_AROUND_TADS", args[1], args[2])
dir.create(outFolder, recursive=TRUE)

inFold <- file.path(paste0("CREATE_SAMPLE_AROUND_TADS_", toupper(args[1])), args[2])
stopifnot(dir.exists(inFold))

corrMeth <- "pearson"
script0_name <- "0_prepGeneData"

pipFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))


all_ds_sample_data_file <- file.path(inFold, paste0("all_ds_sample_around_TADs_", args[1], ".Rdata") )
stopifnot(file.exists(all_ds_sample_data_file))
all_ds_sample_data <- eval(parse(text = load(all_ds_sample_data_file)))

all_ds <- names(all_ds_sample_data)

nDS <- length(all_ds)

cat("... found nDS = ", nDS, "\n")

ds = all_ds[2]
all_ds_around_TADs_corr <- foreach(ds = all_ds ) %do% {
  
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
  
  all_regs <- names(ds_sample_data)
  
  stopifnot(setequal(all_regs, regionList))
  
  reg = all_regs[7]
  ds_all_corr_data <- foreach(reg = all_regs) %dopar% {
    
    tad_data <- ds_sample_data[[paste0(reg)]]
    
    if(tad_data$nGenes > 1) {
      
      sample_genes <- tad_data$genes
      stopifnot(sample_genes %in% geneList)
      
      toKeep <- which(geneList %in% sample_genes)
      
      sub_normDT <- norm_rnaseqDT[toKeep,]  # !!! can be done like this because reordered here above
      
      stopifnot(nrow(sub_normDT) == tad_data$nGenes)
      
      stopifnot(rownames(sub_normDT) == names(geneList)[toKeep])
      stopifnot(cond1_ID %in% colnames(sub_normDT))
      stopifnot(cond2_ID %in% colnames(sub_normDT))
      
      sub_normDT_cond1 <- sub_normDT[,cond1_ID]
      sub_normDT_cond2 <- sub_normDT[,cond2_ID]
      
      stopifnot(nrow(sub_normDT) == nrow(sub_normDT_cond1))
      stopifnot(nrow(sub_normDT) == nrow(sub_normDT_cond2))
      
      stopifnot( ncol(sub_normDT_cond1) + ncol(sub_normDT_cond2) == ncol(sub_normDT))
      stopifnot( ncol(sub_normDT_cond1) == length(cond1_ID))
      stopifnot(ncol(sub_normDT_cond2) == length(cond2_ID))
      
      coexprMat <- cor(t(sub_normDT), method = corrMeth)
      stopifnot(dim(coexprMat) == tad_data$nGenes)
      meanCorr_all <- mean(coexprMat[lower.tri(coexprMat, diag = withDiag)])
      stopifnot(!is.na(meanCorr_all))
      
      coexprMat_cond1 <- cor(t(sub_normDT_cond1), method = corrMeth)
      stopifnot(dim(coexprMat_cond1) == tad_data$nGenes)
      meanCorr_cond1 <- mean(coexprMat_cond1[lower.tri(coexprMat_cond1, diag = withDiag)])
      stopifnot(!is.na(meanCorr_cond1))
      
      coexprMat_cond2 <- cor(t(sub_normDT_cond2), method = corrMeth)
      stopifnot(dim(coexprMat_cond2) == tad_data$nGenes)
      meanCorr_cond2 <- mean(coexprMat_cond2[lower.tri(coexprMat_cond2, diag = withDiag)])
      stopifnot(!is.na(meanCorr_cond2))
      
      
    }else {
      meanCorr_all <- NA
      meanCorr_cond1 <- NA
      meanCorr_cond2 <- NA
    }
    
    list(
      nGenes = tad_data$nGenes,
      meanCorr = meanCorr_all,
      meanCorr_cond1 = meanCorr_cond1,
      meanCorr_cond2 = meanCorr_cond2
    )
    
  } # end iterating over all TADs for the current dataset
  names(ds_all_corr_data) <- all_regs
  ds_all_corr_data
} # end iterating over all DS
names(all_ds_around_TADs_corr) <- all_ds 
  
outFile <- file.path(outFolder, "all_ds_around_TADs_corr.Rdata")
save(all_ds_around_TADs_corr, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

  
  
####################################################################################
####################################################################################3
####################################################################################3
txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))




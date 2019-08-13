options(scipen=100)

setDir <- ""

# Rscript check_TCGA_exprMinCpm_nbr.R

script_name <- "check_TCGA_exprMinCpm_nbr.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(lattice)

registerDoMC(40)

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

script0_name <- "0_prepGeneData"

# % of genes at least 5 reads in x % of samples
mainFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

outFolder <- "CHECK_TCGA_EXPRMINCPM_NBR"
dir.create(outFolder, recursive = TRUE)

hicds = all_hicds[1]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
# check if in pipeline region list
exprds = all_exprds[[paste0(hicds)]]

minReads <- 5

minSample <- 0.1
minSampleRatio_seq <- seq(0, 1, 0.05)


all_geneNbr_dt <- foreach(hicds = all_hicds, .combine='cbind') %dopar% {
  
  exprds_generatio_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='cbind') %do% {
    

  rnadt_file <- file.path(pipFolder, hicds, exprds, script0_name, "rna_rnaseqDT.Rdata")
  stopifnot(file.exists(rnadt_file))
  rna_DT <- eval(parse(text = load(rnadt_file)))
  rnagene_file <- file.path(pipFolder, hicds, exprds, script0_name, "rna_geneList.Rdata")
  stopifnot(file.exists(rnagene_file))
  rna_gene <- eval(parse(text = load(rnagene_file)))
  stopifnot(names(rna_gene) %in% rownames(rna_DT))
  rna_DT <- rna_DT[names(rna_gene),]
  stopifnot(nrow(rna_DT) == length(rna_gene))
  
  totSamples <- ncol(rna_DT)
  
  stopifnot(rna_DT >= 0)
  
  source(file.path(mainFolder, "PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R")))
  samp1 <- eval(parse(text=load(sample1_file)))
  samp2 <- eval(parse(text=load(sample2_file)))
  stopifnot(totSamples == length(samp1) + length(samp2))
  
  genes_nSample_atLeastMinReads <- apply(rna_DT, 1, function(x) sum(x >= minReads))
  stopifnot(length(genes_nSample_atLeastMinReads) == nrow(rna_DT))
  
  nbrGenes_withMinSample <- sapply(minSampleRatio_seq, function(min_sampleNbr_ratio) {
    min_sampleNbr <- min_sampleNbr_ratio * totSamples
    # sum(genes_nSample_atLeastMinReads >= min_sampleNbr)/nrow(rna_DT) # => the ratio
    sum(genes_nSample_atLeastMinReads >= min_sampleNbr) # the actual number
  })
  names(nbrGenes_withMinSample) <- minSampleRatio_seq
  
  # plot(
  #   x = minSampleRatio_seq,
  #   y = nbrGenes_withMinSample,
  #   ylab = paste0("% of genes passing threshold"),
  #   xlab = paste0("threshold % of samples with at least ", minReads, " reads (RSEM)")
  # )
  
  
  dataset <- paste0(hicds, "-", exprds)
  dt <- data.frame(y = nbrGenes_withMinSample)
  colnames(dt) <- dataset
  rownames(dt) <- names(nbrGenes_withMinSample)
  dt
  }
  exprds_generatio_dt
}
outFile <- file.path(outFolder, "all_geneNbr_dt.Rdata")
save(all_geneNbr_dt, file=outFile, version=2)  
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("geneRatio_vs_sampleRatio.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
plot(NULL, 
     xlim=c(0,1), 
     ylim=c(min(all_geneNbr_dt),max(all_geneNbr_dt)),
     ylab = paste0("Nbr genes passing threshold"),
     xlab = paste0("Ratio samples with at least ", minReads, " counts"),# (RSEM)"),
     axes = F
)
xcoord <- as.numeric(rownames(all_geneNbr_dt))
stopifnot(!is.na(xcoord))
foo <- apply(all_geneNbr_dt, 2, function(x) lines(x=xcoord, y = x))
foo <- sapply(seq(0.8, 1, by=0.05), function(x) abline(v=x, lty=2, col="darkgrey"))
axis(2)
axis(1, at = xcoord, labels = xcoord, las=2)
box(bty="l")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


#######################################################################################################################################


v0_nbr_genes <- foreach(hicds = all_hicds, .combine='c') %dopar% {
  
  exprds_genes <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='c') %do% {
    gene_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(gene_file))
    geneList <- eval(parse(text = load(gene_file)))
    
    # rnagene_file <- file.path(pipFolder, hicds, exprds, script0_name, "rna_geneList.Rdata")
    # stopifnot(file.exists(rnagene_file))
    # rna_gene <- eval(parse(text = load(rnagene_file)))
    
    # length(geneList)/length(rna_gene) # the ratio
    length(geneList) # the actual number
  }
  names(exprds_genes) <- paste0(hicds, "-", all_exprds[[paste0(hicds)]])
  exprds_genes
}

thresh <- 0.9

v2_nbr_genes <- setNames(all_geneNbr_dt[paste0(thresh),], colnames(all_geneNbr_dt))

stopifnot(setequal(names(v2_nbr_genes), names(v0_nbr_genes)))

plot(x = v0_nbr_genes,
     y = v2_nbr_genes[names(v0_nbr_genes)])



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

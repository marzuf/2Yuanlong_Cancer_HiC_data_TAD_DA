# Rscript corr_around_BD.R

# for a given number of genes
# go at each BD location
# select "n" genes left and right of the BD
# to get an empirical distribution of the correlation
# between expression of genes separated by the boundary

script_name <- "corr_around_BD.R"

cat("... start ", script_name, "\n")

startTime <- Sys.time()

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


plotType <- "png"
myHeight <- 400
myWidth <- 600

SSHFS <- FALSE
nCpu <- ifelse(SSHFS, 2, 40)

buildCorrAroundBD <- TRUE
buildObsCorr <- TRUE

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(nCpu)

outFolder <- file.path("CORR_AROUND_BD")
dir.create(outFolder)

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

stopifnot(dir.exists(pipOutFolder))

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(length(all_hicexpr_ds) > 0)
stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))

ds="ENCSR862OGI_RPMI-7951_40kb/TCGAskcm_lowInf_highInf"
ds=all_hicexpr_ds[1]

corMet <- "pearson"
withDiago <- FALSE

minGeneNbr <- 3

sideNbrMax <- 10
nMaxGenes <- 2*sideNbrMax

nameVec <- paste0("nGenes", 1:sideNbrMax)

all_gene_nbrs <- 1:sideNbrMax

all_genesTot_wilcox <- paste0("totGenes", (minGeneNbr):(nMaxGenes))



all_ds_corrData <- list()

if(buildCorrAroundBD) {
  # all_hicexpr_ds=all_hicexpr_ds[1]
  for(ds in all_hicexpr_ds) {
    hicds <- file.path(dirname(ds))
    exprds <- basename(ds)
    stopifnot(dir.exists(hicds))
    
    cat("... start dataset: ", exprds, "\n")
    
    dsPipOutDir <- file.path(pipOutFolder, ds)
    stopifnot(dir.exists(dsPipOutDir))
    
    ### RETRIEVE THE GENE2TAD ASSIGNMENT
    g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    
    ### RETRIEVE THE TAD POSITIONS
    tadposFile <- file.path(hicds, "genes2tad", "all_assigned_regions.txt")
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
    
    
    bdpos_DT1 <- tadpos_DT[, c("chromo", "start")]
    bdpos_DT1$start <- bdpos_DT1$start - 1
    colnames(bdpos_DT1) <- c("chromo", "BDpos")
    bdpos_DT2 <- tadpos_DT[, c("chromo", "end")]
    colnames(bdpos_DT2) <- c("chromo", "BDpos")
    bdpos_DT <- rbind(bdpos_DT1, bdpos_DT2)
    stopifnot(nrow(bdpos_DT) == 2*nrow(tadpos_DT))
    bdpos_DT <- unique(bdpos_DT)
    nPos <- nrow(bdpos_DT)
    
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
    
    norm_rnaseqDT <- eval(parse(text = load(file.path(dsPipOutDir, script0_name, "rna_qqnorm_rnaseqDT.Rdata"))))
    stopifnot(names(pipeline_geneList) %in% rownames(norm_rnaseqDT))
    # stopifnot(pipeline_geneList %in% rownames(norm_rnaseqDT)) -> FALSE
    norm_rnaseqDT <- norm_rnaseqDT[names(pipeline_geneList),]    
    stopifnot(rownames(norm_rnaseqDT) == names(pipeline_geneList))
    
    # nPos=10
    allBD_nGenes_corr <- foreach(i_bd = seq_len(nPos)) %dopar% {
    # allBD_nGenes_corr <- foreach(i_bd = seq_len(nPos)) %do% {
      
      cat("...... start BD pos. : \t", i_bd, "/", nPos, "\n")
      
      
      curr_chromo <- as.character(bdpos_DT$chromo[i_bd])
    
      curr_pos <- bdpos_DT$BDpos[i_bd]
      stopifnot(is.numeric(curr_pos))
      stopifnot(length(curr_pos) == 1)
      
      # !!! EXTRACT GENES BASED ON START POSITION RELATIVE TO BD
      # !!! SMALLER THAN / GREATER *OR EQUAL* THAN BD POSITION (smaller not equal otherwise genes could come twice)
      
      curr_g2t <- g2t_DT[g2t_DT$chromo == curr_chromo,,drop=FALSE]
      stopifnot(nrow(curr_g2t) > 0)
      
      stopifnot(is.numeric(curr_g2t$start), is.numeric(curr_g2t$end))
      curr_g2t <- curr_g2t[order(curr_g2t$start, curr_g2t$end),,drop=FALSE]
      
      
      curr_g2t$posDiff <- curr_pos - curr_g2t$start
      
      curr_genesLeftDT <- curr_g2t[curr_g2t$start < curr_pos,,drop=FALSE] 
      stopifnot(nrow(curr_genesLeftDT) == sum(curr_g2t$posDiff > 0))
      
      curr_genesRightDT <- curr_g2t[curr_g2t$start >= curr_pos,,drop=FALSE] 
      stopifnot(nrow(curr_genesRightDT) == sum(curr_g2t$posDiff <= 0))
      
      nrow(curr_genesRightDT)
      nrow(curr_genesLeftDT)
      
      if(nrow(curr_genesRightDT) == 0 | nrow(curr_genesLeftDT) == 0) {
        return(NA)
      }
  
      stopifnot(curr_genesRightDT$entrezID %in% pipeline_geneList)   
      stopifnot(curr_genesLeftDT$entrezID %in% pipeline_geneList)   
      
      # on the right: negative, sort from biggest to smallest (decreasing)
      # on the left: positive, sort from smallest to biggest (increasing)
      curr_genesRightDT <- curr_genesRightDT[order(curr_genesRightDT$posDiff, decreasing = TRUE),,drop=FALSE]
      curr_genesLeftDT <- curr_genesLeftDT[order(curr_genesLeftDT$posDiff, decreasing = FALSE),,drop=FALSE]
      
      stopifnot(is.numeric(curr_genesLeftDT$posDiff))
      stopifnot(is.numeric(curr_genesRightDT$posDiff))
      
      # > all(pipeline_geneList %in% rownames(rna_qqnorm_rnaseqDT))
      # [1] FALSE
      # > all(names(pipeline_geneList) %in% rownames(rna_qqnorm_rnaseqDT))
      # [1] TRUE
      
      # stopifnot(as.character(curr_genesRightDT$entrezID) %in% as.character(rownames(rna_qqnorm_rnaseqDT)))
      # stopifnot(as.character(curr_genesLeftDT$entrezID) %in% as.character(rownames(rna_qqnorm_rnaseqDT)))
      
      stopifnot(names(pipeline_geneList)[pipeline_geneList %in% as.character(curr_genesRightDT$entrezID)] %in%
                  as.character(rownames(rna_qqnorm_rnaseqDT)))
      stopifnot(names(pipeline_geneList)[pipeline_geneList %in% as.character(curr_genesLeftDT$entrezID)] %in% 
                  as.character(rownames(rna_qqnorm_rnaseqDT)))
      
      #all_gene_nbrs=1:10
      
      bd_all_nGenes_corr <- rep(NA, length(all_gene_nbrs))
      names(bd_all_nGenes_corr) <- as.character(all_gene_nbrs)
  
      for(nGenes in all_gene_nbrs) {
        
        if( nrow(curr_genesLeftDT) < nGenes | nrow(curr_genesRightDT) < nGenes ) {
          bd_all_nGenes_corr[as.character(nGenes)] <- NA
          next
        }
        
        cat("......... start # genes = ", nGenes, "\n")
        
        # select nGenes numbers from each side of each boundary
        # > all(curr_g2t$entrezID %in% names(pipeline_geneList))
        # [1] FALSE
        # > all(curr_g2t$entrezID %in% (pipeline_geneList))
        # [1] TRUE
        # 
        # entrezID <-> geneList
        # names(geneList) <-> rownames rnaseq
        
        stopifnot(curr_genesLeftDT$entrezID %in% pipeline_geneList)
        stopifnot(curr_genesRightDT$entrezID %in% pipeline_geneList)
        
        nGenesLeft <- curr_genesLeftDT$entrezID[1:nGenes]
        stopifnot(nGenesLeft %in% pipeline_geneList)
        # stopifnot(nGenesLeft %in% names(pipeline_geneList))
        nGenesLeft_rnaID <- names(pipeline_geneList)[pipeline_geneList %in% nGenesLeft]
        stopifnot(nGenesLeft_rnaID %in% rownames(rna_qqnorm_rnaseqDT))
        
        stopifnot(!is.na(nGenesLeft))
        stopifnot(is.character(nGenesLeft))
        
        nGenesRight <- curr_genesRightDT$entrezID[1:nGenes]
        stopifnot(nGenesRight %in% pipeline_geneList)
        # stopifnot(nGenesRight %in% names(pipeline_geneList))
        nGenesRight_rnaID <- names(pipeline_geneList)[pipeline_geneList %in% nGenesRight]
        stopifnot(nGenesRight_rnaID %in% rownames(rna_qqnorm_rnaseqDT))
        
        stopifnot(!nGenesRight %in% nGenesLeft)
        
        stopifnot(!is.na(nGenesRight))
        stopifnot(is.character(nGenesRight))
        
        stopifnot(length(nGenesLeft) == nGenes)
        stopifnot(length(nGenesRight) == nGenes)
      
        ### OR PUT AS A SINGLE VECTOR ???
        exprLeft <- rna_qqnorm_rnaseqDT[as.character(nGenesLeft_rnaID),,drop=FALSE]
        exprRight <- rna_qqnorm_rnaseqDT[as.character(nGenesRight_rnaID),,drop=FALSE]
        stopifnot(nrow(exprLeft) == nrow(exprRight))
        stopifnot(dim(exprLeft) == dim(exprRight))
        # cor(as.numeric(exprLeft), as.numeric(exprRight))
            
        exprRightLeft <- rna_qqnorm_rnaseqDT[c(as.character(nGenesLeft_rnaID), as.character(nGenesRight_rnaID)),,drop=FALSE]
        stopifnot(nrow(exprRightLeft) == nrow(exprLeft) + nrow(exprRight))
      
        corMatrixRightLeft <- cor(t(exprRightLeft), method=corMet)
        stopifnot(dim(corMatrixRightLeft) == nGenes*2)
        
        meanCorr_rightLeft <- mean(corMatrixRightLeft[lower.tri(corMatrixRightLeft, diag = withDiago)], na.rm=TRUE)
        
        bd_all_nGenes_corr[as.character(nGenes)] <- meanCorr_rightLeft # should match the name used when create vector full of NA
        
      } # end for-iterating over gene #
      stopifnot(length(bd_all_nGenes_corr) == length(all_gene_nbrs))
      #stopifnot(!is.na(bd_all_nGenes_corr)) # not true if not enough genes on the right or on the right
      
      names(bd_all_nGenes_corr) <- paste0("nGenes", names(bd_all_nGenes_corr))
      bd_all_nGenes_corr
      
    } # end foreach-iterating over boundary positions
    stopifnot(length(allBD_nGenes_corr) == nPos)
    names(allBD_nGenes_corr) <- paste0("boundary", 1:nPos)
    all_ds_corrData[[paste0(ds)]] <- allBD_nGenes_corr
  } # end for iterating over ds
  
    
  outFile <- file.path(outFolder, "all_ds_corrData.Rdata")
  save(all_ds_corrData, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))

# end-if buildCorrAroundBD
} else {
  outFile <- file.path(outFolder, "all_ds_corrData.Rdata")
  all_ds_corrData <- eval(parse(text = load(outFile)))
}

across_nDS <- length(all_ds_corrData)

######### INVESTIGATE ACROSS BD DATA

# look at the result:
# the data look like
# $`ENCSR079VIJ_G401_40kb/TCGAkich_norm_kich`$boundary977
# nGenes1     nGenes2     nGenes3     nGenes4     nGenes5     nGenes6     nGenes7     nGenes8 
# -0.04657960 -0.01750222  0.05989178  0.08549617  0.11009050  0.14444184  0.10880600  0.12563351 
# nGenes9    nGenes10 
# 0.13755671  0.11581646 


######### PREPARE OBSERVED DATA
script0_name <- "0_prepGeneData"
script4_name <- "4_runMeanTADCorr"

ds = names(all_ds_corrData)[1]

if(buildObsCorr) {

  obs_corrData <- foreach(ds = names(all_ds_corrData)) %do% {
    
    ds_outFolder <- file.path("PIPELINE", "OUTPUT_FOLDER", ds)
    
    ### LOAD THE CORRELATION DATA FOR THIS DATASET
    obs_meanCorrFile <- file.path(ds_outFolder, script4_name, "all_meanCorr_TAD.Rdata")
    stopifnot(file.exists(obs_meanCorrFile))
    obs_meanCorr <- eval(parse(text=load(obs_meanCorrFile)))
    
    ### KEEP ONLY THE TADs USED IN THE PIPELINE
    stopifnot(dir.exists(file.path(ds_outFolder, script0_name)))
    tadListFile <- file.path(ds_outFolder, script0_name, "pipeline_regionList.Rdata")
    stopifnot(file.exists(tadListFile))
    pipeline_tadList <- eval(parse(text = load(tadListFile))) # not adjusted
    
    ### RETRIEVE THE GENES USED IN THE PIPELINE - script0
    stopifnot(dir.exists(file.path(ds_outFolder, script0_name)))
    geneListFile <- file.path(ds_outFolder, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneListFile))
    pipeline_geneList <- eval(parse(text = load(geneListFile))) 
    
    g2tFile <- file.path(dirname(ds), "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    
    stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
    
    g2t_DT <- g2t_DT[g2t_DT$entrezID %in% pipeline_geneList,]
    
    tadSize <- setNames(as.numeric(table(g2t_DT$region)), names(table(g2t_DT$region)))
    
    stopifnot(pipeline_tadList %in% names(tadSize))
    
    tadSize <- tadSize[names(tadSize) %in% pipeline_tadList]
    stopifnot(min(tadSize) >= minGeneNbr)
    
    all_sizes <- sort(unique(tadSize))
    
    tsize = all_sizes[1]
    corr_by_size <- foreach(tsize = all_sizes) %do% {
      curr_tads <- names(tadSize)[tadSize == tsize]
      stopifnot(length(curr_tads) > 0)
      # retrieve the correlation for these TADs
      stopifnot(curr_tads %in% names(obs_meanCorr))
      obs_meanCorr[names(obs_meanCorr) %in% curr_tads]
    }
    names(corr_by_size) <- paste0("totGenes", all_sizes)
    
    corr_by_size
    
  }
  names(obs_corrData) <- names(all_ds_corrData)
  outFile <- file.path(outFolder, "obs_corrData.Rdata")
  save(obs_corrData, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else{ # end-if buildObsCorr
  outFile <- file.path(outFolder, "obs_corrData.Rdata")
  obs_corrData <- eval(parse(text = load(outFile)))
}

obs_nDS <- length(obs_corrData)

########################################################
######################################################## PLOT INTRA TAD VS ACROSS DATA 
########################################################

cat("... start plotting:\tintraTAD vs. observed data only\n")

all_genesTot <- unique(unlist(lapply(obs_corrData, names)))

geneTot <- "totGenes3"
for(geneTot in all_genesTot) {
  
  nTotGenes <-as.numeric(gsub("totGenes", "", geneTot))
  stopifnot(!is.na(nTotGenes))
  
  if(nTotGenes > nMaxGenes)
    break
  
  nAcrossGenes <- unique(c(floor(nTotGenes/2), ceiling(nTotGenes/2)))
  acrossName <- paste0("nGenes", nAcrossGenes)
  
  obs_values <- as.numeric(unlist(lapply(obs_corrData, function(x) x[paste0(geneTot)])))
  stopifnot(!is.na(obs_values))
  
  across_values <- unlist(lapply(all_ds_corrData, function(ds_data) {
    lapply(ds_data, function(x) x[names(x) %in% acrossName])
  }))
  

  outFile <- file.path(outFolder, paste0("intraTAD_vs_acrossBD_nTotGenes", nTotGenes,".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(
    list(obs_values=obs_values, 
         across_values = across_values),
    plotTit = paste0("nTotGenes = ", nTotGenes, " (", paste0(nAcrossGenes, collapse="-"), " across)"),
    my_xlab = paste0("intra-TAD/across BD mean corr. (", corMet,")")
  )
  mtext(side=3, text = paste0("(obs_nDS=", obs_nDS, "; across_nDS=", across_nDS, ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}

########################################################
######################################################## wilcox pvals
########################################################


allGenes_allDS_pvals <- foreach(geneTot = all_genesTot_wilcox) %do% {
  
  nTotGenes <-as.numeric(gsub("totGenes", "", geneTot))
  stopifnot(!is.na(nTotGenes))
  
  if(nTotGenes > nMaxGenes)
    break
  
  nAcrossGenes <- unique(c(floor(nTotGenes/2), ceiling(nTotGenes/2)))
  acrossName <- paste0("nGenes", nAcrossGenes)
  
  obs_values_ds <- lapply(obs_corrData, function(x) x[paste0(geneTot)])
  
  across_values_ds <- lapply(all_ds_corrData, function(ds_data) {
    lapply(ds_data, function(x) x[names(x) %in% acrossName])
  })
  
  all_ds <- names(obs_values_ds)
  stopifnot(all_ds %in% names(across_values_ds))
  
  myds=all_ds[1]
  myds=all_ds[40]
  ngenes_allDS_pvals <- foreach(myds = all_ds, .combine='c') %do% {
    wilcox_obs <- unlist(obs_values_ds[[myds]])
    wilcox_perm <- unlist(across_values_ds[[myds]])
    if(is.null(wilcox_obs)) return(NA)
    wilcox.test(
      x=wilcox_obs,
      y=wilcox_perm
    )$p.value
  } # end-foreach iterating over datasets
  names(ngenes_allDS_pvals) <- all_ds
  ngenes_allDS_pvals
} # end-foreach iterating over # of genes
names(allGenes_allDS_pvals) <- all_genesTot_wilcox  
outFile <- file.path(outFolder, "allGenes_allDS_pvals.Rdata")
save(allGenes_allDS_pvals, file = outFile)  
  
allGenes_allDS_pvals_log10 <- lapply(allGenes_allDS_pvals, function(x) -log10(x))

stopifnot(unique(unlist(lapply(allGenes_allDS_pvals_log10, length))) == obs_nDS)

outFile <- file.path(outFolder, paste0("boxplot_wilcoxPvals.", plotType))
do.call(plotType, list(outFile, width=myWidth, height=myHeight))
boxplot(allGenes_allDS_pvals_log10, las=3, main="Wilcox p-val [-log10]")
mtext(side=3, text = paste0("(obs_nDS=", obs_nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
  


########################################################
######################################################## PLOT THE ACROSS DATA ONLY
########################################################

cat("... start plotting:\tacross BD data only\n")


dist_across_meanCorr <- foreach(ngene = nameVec) %do%{
  tmp <- unlist(lapply(all_ds_corrData, function(all_data)
    lapply(all_data, function(x) x[ngene])))
  tmp
}

names(dist_across_meanCorr) <- nameVec
outFile <- file.path(outFolder, paste0("acrossBD_dist_genes1-5.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  dist_across_meanCorr[1:5],
  plotTit = paste0("across BD corr. - all data"),
  my_xlab = paste0("across BD mean corr. (", corMet,")")
)
mtext(side=3, text = paste0("(across_nDS=", across_nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("acrossBD_dist_genes6-10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  dist_across_meanCorr[6:10],
  plotTit = paste0("across BD corr. - all data"),
  my_xlab = paste0("across BD mean corr. (", corMet,")")
)
mtext(side=3, text = paste0("(across_nDS=", across_nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("acrossBD_dist_genes11-15.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  dist_across_meanCorr[11:15],
  plotTit = paste0("across BD corr. - all data"),
  my_xlab = paste0("across BD mean corr. (", corMet,")")
)
mtext(side=3, text = paste0("(across_nDS=", across_nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

########################################################
######################################################## PLOT THE OBSERVED DATA ONLY
########################################################

cat("... start plotting:\tobserved data only\n")

nDS <- length(obs_corrData)

geneTot = "totGenes3"
dist_meanCorr <- foreach(geneTot = all_genesTot) %do%{
  # tmp <- lapply(obs_corrData, function(all_data)
  #   all_data[geneTot])
  tmp <- unlist(lapply(obs_corrData, function(all_data)
    all_data[geneTot]))
    #lapply(all_data, function(x) x[ngene])))
  tmp
}
names(dist_meanCorr) <- all_genesTot

outFile <- file.path(outFolder, paste0("obs_intraCorr_dist_genes1-5.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  dist_meanCorr[1:5],
  plotTit = paste0("intra-TAD mean corr. - all data"),
  my_xlab = paste0("intra-TAD mean corr. (", corMet,")")
)
mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("obs_intraCorr_dist_genes6-10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  dist_meanCorr[6:10],
  plotTit = paste0("intra-TAD mean corr. - all data"),
  my_xlab = paste0("intra-TAD mean corr. (", corMet,")")
)
mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("obs_intraCorr_dist_genes11-15.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  dist_meanCorr[11:15],
  plotTit = paste0("intra-TAD mean corr. - all data"),
  my_xlab = paste0("intra-TAD mean corr. (", corMet,")")
)
mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("obs_intraCorr_dist_genes16-20.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  dist_meanCorr[16:20],
  plotTit = paste0("intra-TAD mean corr. - all data"),
  my_xlab = paste0("intra-TAD mean corr. (", corMet,")")
)
mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("obs_intraCorr_dist_genes21-25.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  dist_meanCorr[21:25],
  plotTit = paste0("intra-TAD mean corr. - all data"),
  my_xlab = paste0("intra-TAD mean corr. (", corMet,")")
)
mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("obs_intraCorr_dist_genes26-30.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  dist_meanCorr[26:30],
  plotTit = paste0("intra-TAD mean corr. - all data"),
  my_xlab = paste0("intra-TAD mean corr. (", corMet,")")
)
mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


########################################################
########################################################
########################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))











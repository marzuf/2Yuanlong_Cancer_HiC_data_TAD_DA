#!/usr/bin/Rscript

# CORRESPONDS TO THE 19_SAMP_emp_measurement.R in TAD_DE_pipeline_v2
# adapted here to use different values for the meanCorr permutation

# in 19_SAMP_emp_measurement.R , use get_SAM_FDR function defined in TAD_DE_pipeline_v2/TAD_DE_utils.R

startTime <- Sys.time()

suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# empirical FDR
# logFC as effect size (measure of actual difference I want)

# rank by logFC
# set a cut-off k above which I say as significant
# N = number of real solutions with logFC > k

# for each permutation, compute how many solutions greater than k
# R = average number of random solutions with logFC > k
# => represents the number of false positives
# => FDR = R/N

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

source(paste0(setDir, "/mnt/ed4/marie/scripts/RNA_seq_v2_before0405/RNAseq_fct.R"))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

script0_name <- "0_prepGeneData"


script4_name <- "4_runMeanTADCorr"
script6_name <- "6_runPermutationsMeanLogFC"
script7_name <- "7_runPermutationsMeanTADCorr"
script8_name <- "8c_runAllDown"

script_name <- "SAM_emp_measurement_meanCorr.R"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

# cat(paste0("setDir = ", setDir, "\n"))
source("main_settings.R") # setDir is the main_settings not in run_settings
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))


# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

# caller <- "DI"
# curr_dataset <- "TCGAbrca_lum_bas"
# outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2", "SAM_emp_measurement", caller, curr_dataset)
# system(paste0("mkdir -p ", outFold))

# TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_", caller, "_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    
# # file with assignment from entrez to all regions
# gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_", caller, "_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_genes_positions.txt") 
# curr_outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER", curr_dataset)

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

fixCutOffSeq <- TRUE

#**************************************************************************************************** COMBINED LOG FC AND INTRA TAD CORR

TADpos_DT <- read.delim(TADpos_file, header=F, stringsAsFactors = F, col.names=c("chromo", "region", "start", "end"))
gene2tad_DT <- read.delim(gene2tadDT_file, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))



############# intraCorr
# obs_intraCorr_file <- file.path(curr_outFold, script4_name, "all_meanCorr_TAD.Rdata")
obs_intraCorr_file <- file.path(pipOutFold, script4_name, "all_meanCorr_TAD.Rdata")
stopifnot(file.exists(obs_intraCorr_file))
obs_vect_intraCorr <- eval(parse(text = load(obs_intraCorr_file)))

# shuff_intraCorr_file <- file.path(curr_outFold, script7_name, "meanCorr_permDT.Rdata")
shuff_intraCorr_file <- file.path(pipOutFold, script7_name, "meanCorr_permDT.Rdata")
stopifnot(file.exists(shuff_intraCorr_file))
permutDT_intraCorr <- eval(parse(text = load(shuff_intraCorr_file)))

cut_off_seq_intraCorr <- seq(0,1,0.1)



###########################################################################################





commonReg <- sort(names(obs_vect_intraCorr))

# geneList_file <- file.path(curr_outFold, script0_name, "pipeline_geneList.Rdata")
geneList_file <- file.path(pipOutFold, script0_name, "pipeline_geneList.Rdata")
stopifnot(file.exists(geneList_file))
geneList <- eval(parse(text = load(geneList_file)))

g2t_DT <- gene2tad_DT[gene2tad_DT$entrezID %in% geneList,]
stopifnot(nrow(g2t_DT) > 0)
nbrGenes <- setNames(as.numeric(table(g2t_DT$region)), names(table(g2t_DT$region)))


nbrGenes <- nbrGenes[commonReg]
stopifnot(length(nbrGenes) > 0)

ii <- cut(nbrGenes, breaks = seq(min(nbrGenes), max(nbrGenes), len = 100),
          include.lowest = TRUE)
colors <- colorRampPalette(c("blue","white", "red"))(99)[ii]




permutDT_intraCorr <- permutDT_intraCorr[commonReg,]
permutDT_ratioDown <- permutDT_ratioDown[commonReg,]
permutDT_prodSignedRatio <- permutDT_prodSignedRatio[commonReg,]




empFDR_list <- list()






#**************************************************************************************************** EMP. FDR FROM INTRA TAD CORRELATION

curr_variable <- "intraTADcorr"
obs_vect <- obs_vect_intraCorr
permutDT <- permutDT_intraCorr
curr_variable_plotName <- curr_variable

interRegion <- intersect(names(obs_vect), rownames(permutDT))
stopifnot(setequal(interRegion, commonReg))
obs_vect <- obs_vect[interRegion]
permutDT <- permutDT[interRegion,]
# cut_off_seq <- round(seq(0, max(obs_vect[-which.max(obs_vect)]), length.out = 10),2)
# cut_off_seq <- seq(0, 0.5, by =0.05)
if(fixCutOffSeq) {
  cut_off_seq <- cut_off_seq_intraCorr
} else{
  cut_off_seq <- round(seq(0, max(obs_vect[-which.max(obs_vect)]), length.out = 10),2)
}
# higher: sum(obs_vect >= cut_off)
empFDR_seq <- unlist(sapply(cut_off_seq, function(x) get_SAM_FDR(obs_vect, permutDT, cut_off = x, symDir = "higher", withPlot = F)))

toKeep <- !is.infinite(empFDR_seq) & !is.na(empFDR_seq)
slopeFDR <- as.numeric(coef(lm(empFDR_seq[toKeep] ~ cut_off_seq[toKeep]))["cut_off_seq[toKeep]"])

# ! higher
nbrObservedSignif <- unlist(sapply(cut_off_seq, function(x) sum(obs_vect >= x)))

# TRY VARIABLE CUT-OFFS: plot FDR and nbrObservSignif ~ cut_off
outFile <- paste0(curr_outFold, "/", "FDR_var_cut_off_", curr_variable, ".", plotType)
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_FDR_with_observedSignif(yaxis_empFDR_vect= empFDR_seq,xaxis_cutoff= cut_off_seq, y2_obsSignif= nbrObservedSignif, variableName=curr_variable_plotName, feature_name="TADs")
cat(paste0("... written: ", outFile, "\n"))
foo <- dev.off()

# PLOT ALL VALUES AND AREA PERMUT WITH FDR for the given cut-off
k_cut_off <- 0.2
outFile <- paste0(curr_outFold, "/", "FDR_", k_cut_off, "_cut_off_", curr_variable, ".", plotType)
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
get_SAM_FDR(obs_vect, permutDT, cut_off = k_cut_off, variableName = curr_variable, symDir = "higher", withPlot = TRUE)
# get_SAM_FDR(obs_vect, permutDT, cut_off = k_cut_off, variableName = curr_variable, symDir = "higher", withPlot = TRUE, minQuant=0, maxQuant = 1)
textTAD <- names(obs_vect)[which(obs_vect >= k_cut_off  ) ] # higher !
if(length(textTAD) > 0)
  text(y=obs_vect[textTAD], x = which(names(obs_vect) %in% textTAD), labels = textTAD, pos=2, col="gray")

cat(paste0("... written: ", outFile, "\n"))
foo <- dev.off()

empFDR_list[[paste0("empFDR_", curr_variable)]] <- setNames(empFDR_seq, paste0(cut_off_seq))
empFDR_list[[paste0("nbrSignif_", curr_variable)]] <- setNames(nbrObservedSignif, paste0(cut_off_seq))
empFDR_list[[paste0("slopeEmpFDR_", curr_variable)]] <- slopeFDR

rm(obs_vect)
rm(permutDT)




rm(obs_vect)
rm(permutDT)
##############################
outFile <- paste0(curr_outFold, "/", "empFDR_list.Rdata")
save(empFDR_list, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

tmp <- unlist(empFDR_list[grep("empFDR", names(empFDR_list))])
tmp <- tmp[!is.na(tmp) & !is.infinite(tmp)]
if(any(tmp > 1)) {
  warning("!!! emp. FDR > 1 found !!!\n")
}

##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

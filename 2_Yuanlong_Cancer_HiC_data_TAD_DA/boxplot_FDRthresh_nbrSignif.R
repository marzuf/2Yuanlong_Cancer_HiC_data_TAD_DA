# ... written: MEANCORR_MEANFC_FDR/GSE58752_liver_40kb/TCGAlihc_wt_mutCTNNB1/signifTADs_FC_and_CorrV0.Rdata
# 


pipOutFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "OUTPUT_FOLDER")
inFolder <- file.path("MEANCORR_MEANFC_FDR")

if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}

outFolder <- "MEANCORR_MEANFC_FDR"


hicds = all_hicds[1]
all_dt <- foreach(hicds = all_hicds, .combine='rbind') %do% {
  
  
  
  exprds = all_exprds[[paste0(hicds)]][1]
  all_dt_hicds <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    
    
    curr_inFolder <- file.path(inFolder, hicds, exprds)
    
    
    v0_file <- file.path(curr_inFolder, "signifTADs_FC_and_CorrV0.Rdata")
    stopifnot(file.exists(v0_file))
    
    v0_FDR <- eval(parse(text = load(v0_file)))
    
    v0_DT <- do.call(rbind,lapply(v0_FDR, data.frame))
    v0_DT$FDR_cut_off <- rownames(v0_DT)
    colnames(v0_DT) <- paste0(colnames(v0_DT), "_v0")
    
    
    allcorr_file <- file.path(curr_inFolder, "signifTADs_FC_and_allCorrs.Rdata")
    stopifnot(file.exists(allcorr_file))
    
    allCorrs_FDR <- eval(parse(text = load(allcorr_file)))
    
    
    
    dt_list_allCorrs <- lapply(allCorrs_FDR, function(sublist1) lapply(sublist1, function(sublist2) do.call(rbind, lapply(sublist2,  data.frame))))
    
    
    
  }
  
  all_dt_hicds
  
}


fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
stopifnot(file.exists(fc_file))
stopifnot(file.exists(corr_file))  

tad_fc <- eval(parse(text = load(fc_file)))
tad_corr <- eval(parse(text = load(corr_file)))

all_regs <- names(tad_fc)
stopifnot(setequal(all_regs, names(tad_corr)))

stopifnot(length(tad_fc) == nTADs)
stopifnot(length(tad_corr) == nTADs)

plotDT <- data.frame(
  regions = all_regs,
  meanFC = as.numeric(tad_fc[all_regs]),
  meanCorr = as.numeric(tad_corr[all_regs]),
  stringsAsFactors = FALSE
)


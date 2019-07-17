# compare the FDR~corr thresholds and FDR~FC thresholds

options(scipen=100)

# Rscript meanCorr_meanFC_FDR.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript meanCorr_meanFC_FDR.R   # to run all datasets in one shot

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "meanCorr_meanFC_FDR.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

source("utils_fct.R")


script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script19_name <- "19_SAM_emp_measurement"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4

pipOutFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))


samInFolder <- file.path("SAM_EMP_MEASUREMENT_MEANCORR")
stopifnot(dir.exists(samInFolder))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2 | length(args) == 0)
hicds <- args[1]
exprds <- args[2]


if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}

outFolder <- "MEANCORR_MEANFC_FDR"
dir.create(outFolder, recursive=TRUE)

all_sampTypes <- c("sameNbr", "sameKb", "fixKb")



fc_col <- "dodgerblue1"
corr_col <- "darkolivegreen4"
fdr_thresh <- 0.2

fdr_thresh_seq <- seq(0.05, 0.5, by=0.05)


for(hicds in all_hicds) {
  
  for(exprds in all_exprds[[paste0(hicds)]]) {
    
    
    tad_file <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
    tadList <- eval(parse(text = load(tad_file)))
    
    # ../Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/19_SAM_emp_measurement/empFDR_list.Rdata
    meanFC_FDR_file <- file.path(pipOutFolder, hicds, exprds, script19_name, "empFDR_list.Rdata")
    stopifnot(file.exists(meanFC_FDR_file))
    all_FDR <- eval(parse(text = load(meanFC_FDR_file)))
    # names(all_FDR)
    # [1] "empFDR_logFC"                "nbrSignif_logFC"             "slopeEmpFDR_logFC"           "empFDR_intraTADcorr"        
    # [5] "nbrSignif_intraTADcorr"      "slopeEmpFDR_intraTADcorr"    "empFDR_ratioDown"            "nbrSignif_ratioDown"        
    # [9] "slopeEmpFDR_ratioDown"       "empFDR_prodSignedRatio"      "nbrSignif_prodSignedRatio"   "slopeEmpFDR_prodSignedRatio"
    nTADs <- length(tadList) 
    stopifnot(nTADs == all_FDR[["nbrSignif_logFC"]][[paste0(0)]])
    
    
    fc_fdr <- all_FDR[["empFDR_logFC"]]
    # the 1st FC value for which FDR is smaller than the thresh
    fc_fdr_cut_off <- min(as.numeric(as.character(na.omit(names(fc_fdr)[fc_fdr <= fdr_thresh])))) # will return Inf is none below the thresh
    
    corr_fdr <- all_FDR[["empFDR_intraTADcorr"]]
    corr_fdr_cut_off <- min(as.numeric(as.character(na.omit(names(corr_fdr)[corr_fdr <= fdr_thresh])))) # will return Inf is none below the thresh
    
    
    outFile <- file.path(outFolder, hicds, exprds, paste0("meanCorrV0", "_and_meanFC_FDR_nbrSignif_FDRthresh", fdr_thresh , ".", plotType))
    dir.create(dirname(outFile), recursive = TRUE)
    do.call(plotType, list(outFile, height = myHeight, width = myWidth*1.4))
    
    plot_two_FDR_with_observedSignif(first_var_FDR= fc_fdr,
                                       first_var_nbrSignif = all_FDR[["nbrSignif_logFC"]],
                                                 scd_var_FDR = corr_fdr,
                                     scd_var_nbrSignif = all_FDR[["nbrSignif_intraTADcorr"]],
                                                 first_var_name = "meanFC",
                                                 scd_var_name = "meanCorrV0",
                                     first_var_col = fc_col,
                                     scd_var_col = corr_col,
                                     ya_cut_off = fdr_thresh)
    
    title(main = paste0(hicds, " - ", exprds))
    mtext(text = paste0("nTADs=", nTADs), side = 3, font=3, line=-1)
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
    
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
    
   
    
    outFileSuffix <- file.path(outFolder, hicds, exprds, paste0("meanCorrV0","_vs_", "meanFC _FDRthresh", fdr_thresh))
    plot_meanFC_meanCorr_FDRthresh(
      dataDT=plotDT,
      var1 = "meanFC",
      var2 = "meanCorr",
      annotCol = "regions",
      var2_name = paste0("meanCorrV0"),
      var1_cutoff = fc_fdr_cut_off,
      var2_cutoff = corr_fdr_cut_off,
      abs_var1 = TRUE,
      plotTit = paste0(hicds, " - ", exprds),
      plotSub = paste0("nTADs=", nTADs),
      fileSuffix = outFileSuffix,
      plotType= plotType
    )
    

    
    signifTADs_FC_and_CorrV0 <- lapply(fdr_thresh_seq , function(cutoff_fdr) {
      fc_cut_off <- min(as.numeric(as.character(na.omit(names(fc_fdr)[fc_fdr <= cutoff_fdr]))))
      corr_cut_off <- min(as.numeric(as.character(na.omit(names(corr_fdr)[corr_fdr <= cutoff_fdr]))))
      nSignif <- sum(abs(plotDT$meanFC) >= fc_cut_off & plotDT$meanCorr >= corr_cut_off)
      list(
        fc_cut_off = fc_cut_off,
        corr_cut_off = corr_cut_off,
        nTADs_signif = nSignif
      )
    })
    names(signifTADs_FC_and_CorrV0) <- paste0(fdr_thresh_seq)
    
    outFile <- file.path(outFolder, hicds, exprds, paste0("signifTADs_FC_and_CorrV0.Rdata"))
    save(signifTADs_FC_and_CorrV0, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    samp_type = "sameNbr"    
    signifTADs_bySampType <- foreach(samp_type = all_sampTypes) %dopar% {
      fixSizeKb <- ifelse(samp_type == "fixKb", 1000000, "")
      
      # SAM_EMP_MEASUREMENT_MEANCORR/fixKb/1000000/Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/all_empFDR_seq.Rdata
      
      meanCorr_FDR_file <- file.path(samInFolder, samp_type, fixSizeKb, hicds, exprds, "all_empFDR.Rdata")
      stopifnot(file.exists(meanCorr_FDR_file))
      
      all_meanCorrSample_FDR <- eval(parse(text = load(meanCorr_FDR_file)))
      
      # [1] "sample_meanCorr_allDS"           "sample_meanCorrLeft_allDS"       "sample_meanCorrRight_allDS"     
      # [4] "sample_meanCorrLeftRight_allDS"  "sample_meanCorr_onlyDS"          "sample_meanCorrLeft_onlyDS"     
      # [7] "sample_meanCorrRight_onlyDS"     "sample_meanCorrLeftRight_onlyDS"
      all_corr_types <- names(all_meanCorrSample_FDR)
      corr_type <- "sample_meanCorr_allDS"
      corr_type <- "sample_meanCorr_onlyDS"
     
      
      signifTADs_bySampType_byCorrType <- foreach(corr_type = all_corr_types) %do% {
        
        
        
        corr_fdr <- all_meanCorrSample_FDR[[paste0(corr_type)]][["empFDR"]]
        corr_fdr_cut_off <- min(as.numeric(as.character(na.omit(names(corr_fdr)[corr_fdr <= fdr_thresh])))) # will return Inf is none below the thresh
        
        
        outFile <- file.path(outFolder, hicds, exprds, paste0(samp_type,"_", corr_type,"_and_meanFC_FDR_nbrSignif_FDRthresh", fdr_thresh , ".", plotType))
        dir.create(dirname(outFile), recursive = TRUE)
        do.call(plotType, list(outFile, height = myHeight, width = myWidth*1.4))
        
        plot_two_FDR_with_observedSignif(first_var_FDR= all_FDR[["empFDR_logFC"]],
                                         first_var_nbrSignif = all_FDR[["nbrSignif_logFC"]],
                                         scd_var_FDR =corr_fdr,
                                         scd_var_nbrSignif = all_meanCorrSample_FDR[[paste0(corr_type)]][["nbrSignif"]],
                                         first_var_name = "meanFC",
                                         scd_var_name = gsub("sample_", "",  corr_type),
                                         first_var_col = fc_col,
                                         scd_var_col = corr_col,
                                         ya_cut_off = fdr_thresh)
        title(main = paste0(hicds, " - ", exprds))
        # mtext(text = paste0(gsub("sample_", "", corr_type), " - nTADs=", nTADs), side = 3, font=3, line=-1)
        # mtext(text = paste0("nTADs=", nTADs), side = 3, font=3, line=-1)
        mtext(text = paste0(samp_type, " - nTADs=", nTADs), side = 3, font=3, line=-1)
        foo <- dev.off()
        cat(paste0("... written: ", outFile, "\n"))
        
        
        
        outFileSuffix <- file.path(outFolder, hicds, exprds, paste0(samp_type,"_", corr_type,"_vs_", "meanFC_FDRthresh", fdr_thresh))
        plot_meanFC_meanCorr_FDRthresh(
          dataDT=plotDT,
          var1 = "meanFC",
          var2 = "meanCorr",
          annotCol = "regions",
          var2_name = gsub("sample_", "", corr_type),
          var1_cutoff = fc_fdr_cut_off,
          var2_cutoff = corr_fdr_cut_off,
          abs_var1 = TRUE,
          plotTit = paste0(hicds, " - ", exprds),
          plotSub = paste0(samp_type, " - ", gsub("sample_", "", corr_type), " - nTADs=", nTADs),
          # plotSub = paste0("nTADs=", nTADs),
          fileSuffix = outFileSuffix,
          plotType= plotType,
          legSuppTxt = paste0("FDR thresh.: ",fdr_thresh)
        )
        

        ### MISS THE MEAN CORR TYPE !!!
        
        signifTADs <- lapply(fdr_thresh_seq , function(cutoff_fdr) {
          fc_cut_off <- min(as.numeric(as.character(na.omit(names(fc_fdr)[fc_fdr <= cutoff_fdr]))))
          corr_cut_off <- min(as.numeric(as.character(na.omit(names(corr_fdr)[corr_fdr <= cutoff_fdr]))))
          nSignif <- sum(abs(plotDT$meanFC) >= fc_cut_off & plotDT$meanCorr >= corr_cut_off)
          list(
            fc_cut_off = fc_cut_off,
            corr_cut_off = corr_cut_off,
            nTADs_signif = nSignif
          )
        })
        names(signifTADs) <- paste0(fdr_thresh_seq)
        

        signifTADs
        
        
      } # end iterating the different meanCorr for the current samp_type
      
      names(signifTADs_bySampType_byCorrType) <- all_corr_types
      signifTADs_bySampType_byCorrType      
      
            
    } # end iterating over samptypes
    names(signifTADs_bySampType) <- all_sampTypes
    signifTADs_FC_and_allCorrs <- signifTADs_bySampType
    outFile <- file.path(outFolder, hicds, exprds, paste0("signifTADs_FC_and_allCorrs.Rdata"))
    save(signifTADs_FC_and_allCorrs, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
  } # end iterating over exprds
  
} # end iterating over hicds







##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

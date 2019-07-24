### CREATE TABLE
# -----------------------
#  TAD_ID | genes | meanFC | meanCorr | adj. comb. p-value | signif. at FDR 0.1 ?  | signif. at FDR 0.2  

# Rscript create_final_table.R

options(scipen = 100)

script_name <- "create_final_table.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script19_name <- "19_SAM_emp_measurement"

require(foreach)
require(doMC)
registerDoMC(40)

pipOutFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "OUTPUT_FOLDER")
inFolder <- file.path("MEANCORR_MEANFC_FDR")

outFolder <- "CREATE_FINAL_TABLE"
dir.create(outFolder, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
hicds <- args[1]
exprds <- args[2]

fdr_thresh1 <- 0.1
fdr_thresh2 <- 0.2
corr_type <- "sample_meanCorr_allDS"


samp_type <- "sameNbr"
fixKbSize <- ""

corrEmpFDR_folder <- file.path("SAM_EMP_MEASUREMENT_MEANCORR", samp_type, fixKbSize)
stopifnot(dir.exists(corrEmpFDR_folder))

notAdjCombPval_folder <- file.path("CREATE_EMPPVALCOMB", samp_type, fixKbSize)
stopifnot(dir.exists(notAdjCombPval_folder))
                                

if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}


hicds = all_hicds[1]
all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %do% {
  exprds = all_exprds[[paste0(hicds)]][1]
  hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    cat("... start building DT for :", hicds, " - ", exprds, "\n")
    
    corrEmpFDR_file <- file.path(corrEmpFDR_folder, hicds, exprds, "all_empFDR.Rdata")
    stopifnot(file.exists(corrEmpFDR_file))
    corrEmpFDR <- eval(parse(text = load(corrEmpFDR_file)))
      
    fcEmpFDR_file <- file.path(pipOutFolder, hicds, exprds, script19_name, "empFDR_list.Rdata")
    stopifnot(file.exists(fcEmpFDR_file))
    fc_empFDR <- eval(parse(text = load(fcEmpFDR_file)))
    
    fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
    stopifnot(file.exists(fc_file))
    tad_fc <- eval(parse(text = load(fc_file)))
    all_regs <- names(tad_fc)
    
    
    corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
    stopifnot(file.exists(corr_file))  
    tad_corr <- eval(parse(text = load(corr_file)))
    
    stopifnot(setequal(all_regs, names(tad_corr)))
    
    notAdjCombPval_file <- file.path(notAdjCombPval_folder, hicds, exprds, "meanCorr_meanLogFC_notAdjCombEmpPval.Rdata")
    stopifnot(file.exists(notAdjCombPval_file))
    tad_notAdjCombPval <- eval(parse(text = load(notAdjCombPval_file)))
    stopifnot(setequal(all_regs, names(tad_notAdjCombPval)))
    tad_adjCombPval <- p.adjust(tad_notAdjCombPval, method="BH")
    
    
    
    fc_FDR <- fc_empFDR[["empFDR_logFC"]]
    # the 1st FC value for which FDR is smaller than the thresh
    fc_FDR_co1 <- min(as.numeric(as.character(na.omit(names(fc_FDR)[fc_FDR <= fdr_thresh1])))) 
    # will return Inf is none below the thresh
    fc_FDR_co2 <- min(as.numeric(as.character(na.omit(names(fc_FDR)[fc_FDR <= fdr_thresh2])))) 
    
    corr_FDR <- corrEmpFDR[[paste0(corr_type)]][["empFDR"]]
    corr_FDR_co1 <- min(as.numeric(as.character(na.omit(names(corr_FDR)[corr_FDR <= fdr_thresh1])))) # will return Inf is none below the thresh
    corr_FDR_co2 <- min(as.numeric(as.character(na.omit(names(corr_FDR)[corr_FDR <= fdr_thresh2])))) # will return Inf is none below the thresh
    
    FDR_signif_co1 <- setNames(abs(tad_fc[all_regs]) >= fc_FDR_co1 &
                                 tad_corr[all_regs] >= corr_FDR_co1, all_regs)
    
    stopifnot(setequal(all_regs, names(FDR_signif_co1)))
    
    
    FDR_signif_co2 <- setNames(abs(tad_fc[all_regs]) >= fc_FDR_co2 &
                                 tad_corr[all_regs] >= corr_FDR_co2, all_regs)
    
    stopifnot(setequal(all_regs, names(FDR_signif_co2)))
    
    
    exprds_hicds_dt <- data.frame(
      hicds = hicds,
      exprds = exprds,
      region = all_regs,
      meanLogFC = tad_fc[all_regs],
      meanCorr = tad_corr[all_regs],
      meanLogFC_FDR_co1 = fc_FDR_co1,
      meanLogFC_FDR_co2 = fc_FDR_co2,
      adjPvalComb = tad_adjCombPval[all_regs],
      signifFDRco1 = FDR_signif_co1[all_regs],
      signifFDRco2 = FDR_signif_co2[all_regs],
      stringsAsFactors = FALSE
    )
    rownames(exprds_hicds_dt) <- NULL
    colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "signifFDRco1"] <- paste0("signifFDR_", 
                                                                                     fdr_thresh1)
    colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "signifFDRco2"] <- paste0("signifFDR_", 
                                                                                     fdr_thresh2)
    
    colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "meanLogFC_FDR_co1"] <- paste0("meanLogFC_thresh_FDR", 
                                                                                     fdr_thresh1)
    colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "meanLogFC_FDR_co2"] <- paste0("meanLogFC_thresh_FDR", 
                                                                                     fdr_thresh2)
    
    
    exprds_hicds_dt <- exprds_hicds_dt[order(exprds_hicds_dt$adjPvalComb),]
    
    exprds_hicds_dt
  } # end-foreach iterating over exprds
  hicds_dt
} # end-foreach iterating over hicds
outFile <- file.path(outFolder, "all_result_dt.Rdata")
save(all_result_dt, file = outFile)
cat(paste0("... written: ", outFile,  "\n"))
    

#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))


    
    
      
      
    

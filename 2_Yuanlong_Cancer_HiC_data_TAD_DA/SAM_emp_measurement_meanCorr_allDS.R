
# Rscript SAM_emp_measurement_meanCorr_allDS.R 

# CORRESPONDS TO THE 19_SAMP_emp_measurement.R in TAD_DE_pipeline_v2
# adapted here to use different values for the meanCorr permutation

# in 19_SAMP_emp_measurement.R , use get_SAM_FDR function defined in TAD_DE_pipeline_v2/TAD_DE_utils.R

startTime <- Sys.time()

options(scipen=100)


# empirical FDR

# set a cut-off k above which I say as significant
# N = number of real solutions with logFC > k

# for each permutation, compute how many solutions greater than k
# R = average number of random solutions with logFC > k
# => represents the number of false positives
# => FDR = R/N

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

require(foreach)
source("utils_fct.R")

script_name <- "SAM_emp_measurement_meanCorr_allDS.R"


outFold <- "SAM_EMP_MEASUREMENT_MEANCORR_ALLDS"
dir.create(outFold)

all_sampTypes <- c("sameKb", "fixKb", "sameNbr")

pipFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")

# script0_name <- "0_prepGeneData"
script4_name <- "4_runMeanTADCorr"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight * 1.2

k_cut_off <- 0.5

sort_plotQuantileCutOff <- TRUE


all_obs_corr_files <- list.files(pipOutFolder, pattern = "all_meanCorr_TAD.Rdata", full.names = TRUE, recursive = TRUE)
stopifnot(file.exists(all_obs_corr_files))
stopifnot(length(all_obs_corr_files) > 0)
all_obs_corr_values <- foreach(corr_file = all_obs_corr_files, .combine='c') %dopar% {
  hicds <- basename(dirname(dirname(dirname(corr_file))))
  exprds <- basename(dirname(dirname(corr_file)))
  corr_values <- eval(parse(text = load(corr_file)))
  names(corr_values) <- paste0(hicds, "-", exprds, "\n", names(corr_values))
  corr_values
}

samp_type="sameNbr"
# all_sampTypes=all_sampTypes[1]

foo <- foreach(samp_type = all_sampTypes) %do% {
  
  cat("> START ", samp_type, "\n")
  
    if(samp_type == "fixKb") {
      fixKbSize <- 1000000
    } else{
      fixKbSize <- ""
    }

    currOutFold <- file.path(outFold, samp_type, fixKbSize)
  
    samp_file <- file.path("COEXPR_AROUND_TADS", samp_type, 
                           fixKbSize, "all_ds_around_TADs_corr.Rdata")
    
    stopifnot(file.exists(samp_file))
    
    all_ds_sampleCorr_around <- eval(parse(text = load(samp_file)))
    nbr_permut <- length(all_ds_sampleCorr_around)
    
    sample_meanCorr_allDS <- extract_corr_values(all_ds_sampleCorr_around, "meanCorr", names(all_ds_sampleCorr_around))
    stopifnot(length(sample_meanCorr_allDS) > 0)
    
    sample_meanCorrLeft_allDS <- extract_corr_values(all_ds_sampleCorr_around, "meanCorr_left", names(all_ds_sampleCorr_around))
    stopifnot(length(sample_meanCorrLeft_allDS) > 0)
    
    sample_meanCorrRight_allDS <- extract_corr_values(all_ds_sampleCorr_around, "meanCorr_right", names(all_ds_sampleCorr_around))
    stopifnot(length(sample_meanCorrRight_allDS) > 0)
    
    sample_meanCorrLeftRight_allDS <- c(sample_meanCorrLeft_allDS, sample_meanCorrRight_allDS)
    
    cut_off_seq_intraCorr <- seq(0,1,0.1)
    
    #****************************** EMP. FDR FROM INTRA TAD CORRELATION
    
    all_samplings <- c("sample_meanCorr_allDS","sample_meanCorrLeft_allDS","sample_meanCorrRight_allDS", "sample_meanCorrLeftRight_allDS")
    
    # all_samplings=all_samplings[1]
    all_empFDR_seq <- foreach(sampling_vals_type = all_samplings) %do% {
      
      if(grepl("onlyDS", sampling_vals_type)){
        nbrPermut <- 1
      }else{
        nbrPermut <- nbr_permut 
      }
      
      if(grepl("LeftRight", sampling_vals_type)){
        nbrPermut <- nbrPermut * 2
      } 
      
      
      
      
      cat("...... start: ", sampling_vals_type, "\n")
      
      sampling_values <- eval(parse(text = sampling_vals_type))
      
      cat("...... ", "get_SAM_FDR_aroundTADs", "\n")
      # higher: sum(obs_vect >= cut_off)
      empFDR_seq <- sapply(cut_off_seq_intraCorr, function(x) 
        get_SAM_FDR_aroundTADs(obs_vect = all_obs_corr_values, 
                               permut_values = sampling_values, 
                               cut_off = x, symDir = "higher", 
                               nPermut = nbrPermut,
                               withPlot = F))
      names(empFDR_seq) <- as.character(cut_off_seq_intraCorr)
      
      curr_variable <- sampling_vals_type #"meanTADCorr"
      curr_variable_plotName <- curr_variable
      
      # toKeep <- !is.infinite(empFDR_seq) & !is.na(empFDR_seq)
      # slopeFDR <- as.numeric(coef(lm(empFDR_seq[toKeep] ~ cut_off_seq_intraCorr[toKeep]))["cut_off_seq_intraCorr[toKeep]"])
      
      # ! higher
      nbrObservedSignif <- unlist(sapply(cut_off_seq_intraCorr, function(x) sum(all_obs_corr_values >= x)))
      
      # TRY VARIABLE CUT-OFFS: plot FDR and nbrObservSignif ~ cut_off
      cat("...... ", "plot_FDR_with_observedSignif", "\n")
      outFile <- file.path(currOutFold, paste0("FDR_var_cut_off_", curr_variable, "_", samp_type, ".", plotType))
      dir.create(dirname(outFile), recursive = T)
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      plot_FDR_with_observedSignif(yaxis_empFDR_vect= empFDR_seq,
                                   xaxis_cutoff= cut_off_seq_intraCorr, 
                                   y2_obsSignif= nbrObservedSignif, 
                                   variableName=curr_variable_plotName, 
                                   feature_name="TADs")
      cat(paste0("... written: ", outFile, "\n"))
      foo <- dev.off()
      
      # PLOT ALL VALUES AND AREA PERMUT WITH FDR for the given cut-off
      cat("...... ", "get_SAM_FDR_aroundTADs", "\n")
      outFile <- file.path(currOutFold, paste0("FDR_", k_cut_off, "_cut_off_", curr_variable, "_", samp_type, ".", plotType))
      dir.create(dirname(outFile), recursive = T)
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      get_SAM_FDR_aroundTADs(obs_vect = all_obs_corr_values, 
                             permut_values = sampling_values, 
                             cut_off = k_cut_off, 
                             variableName = curr_variable, 
                             symDir = "higher", 
                             withPlot = TRUE,
                             sortToPlot = sort_plotQuantileCutOff)
      # get_SAM_FDR(obs_vect, permutDT, cut_off = k_cut_off, variableName = curr_variable, symDir = "higher", withPlot = TRUE, minQuant=0, maxQuant = 1)
      textTAD <- names(all_obs_corr_values)[which(all_obs_corr_values >= k_cut_off  ) ] # higher !
      if(length(textTAD) > 0)
        text(y=all_obs_corr_values[textTAD], x = which(names(all_obs_corr_values) %in% textTAD), labels = textTAD, pos=2, col="gray")
      
      cat(paste0("... written: ", outFile, "\n"))
      foo <- dev.off()
      
      empFDR_seq
    }
    names(all_empFDR_seq) <- all_samplings
    
    ##############################
    outFile <- file.path(currOutFold,  "all_empFDR_seq.Rdata")
    dir.create(dirname(outFile))
    save(all_empFDR_seq, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
} # end foreach iterating over sampTypes


##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

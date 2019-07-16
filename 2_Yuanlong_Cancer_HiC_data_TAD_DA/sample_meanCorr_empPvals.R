
# Rscript sample_meanCorr_empPvals.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS
# Rscript sample_meanCorr_empPvals.R 

options(scipen=100)

script_name <- "sample_meanCorr_empPvals.R"
cat("... start ", script_name, "\n")
startTime <- Sys.time()


require(foreach)
require(doMC)

registerDoMC(40)

outFold <- "SAMPLE_MEANCORR_EMPPVALS"
dir.create(outFold)

args <- commandArgs(trailingOnly = TRUE)

hicds <- args[1]
exprds <- args[2]

pipMainFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "OUTPUT_FOLDER")

if(length(args) == 0) {
  all_hicds <- list.files(pipMainFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipMainFolder, x)))
} else {
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}


pipMainFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "OUTPUT_FOLDER")

script4_name <- "4_runMeanTADCorr"


all_sampTypes <- c("sameKb", "fixKb", "sameNbr")
all_corrTypes <- c("meanCorr", "meanCorr_right", "meanCorr_left", "meanCorrLeftRight")


for(hicds in all_hicds) {
  
  for(exprds in all_exprds[[paste0(hicds)]]) {
    
    cat("... start ", hicds, " - ", exprds, "\n")
    
    obs_corr_file <- file.path(pipMainFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
    stopifnot(file.exists(obs_corr_file))
    
    all_obs_corr <- eval(parse(text = load(obs_corr_file)))
    

    samp_type = "sameKb"
    all_cors_empPval_dt <- foreach(samp_type = all_sampTypes, .combine='cbind') %do% {
      
      if(samp_type == "fixKb") {
        fixKbSize <- 1000000
      } else{
        fixKbSize <- ""
      }
      samp_file <- file.path("COEXPR_AROUND_TADS", samp_type, fixKbSize, "all_ds_around_TADs_corr.Rdata")
      stopifnot(file.exists(samp_file))
      all_ds_corr_around <- eval(parse(text = load(samp_file)))
      
      
      corr_type = "meanCorr"
      corr_empPval_dt <- foreach(corr_type = all_corrTypes, .combine='cbind') %dopar% {
        
        
    
        if(corr_type == "meanCorrLeftRight"){
          all_samp_corrs_left <-  unlist(lapply(all_ds_corr_around, function(sub_data){
            lapply(sub_data, function(x) x[[paste0("meanCorr_left")]])
          }))
          all_samp_corrs_left <- na.omit(all_samp_corrs_left)
          all_samp_corrs_right <-  unlist(lapply(all_ds_corr_around, function(sub_data){
            lapply(sub_data, function(x) x[[paste0("meanCorr_right")]])
          }))
          all_samp_corrs_right <- na.omit(all_samp_corrs_left)
          all_samp_corrs <- c(all_samp_corrs_left, all_samp_corrs_right)
        } else {
          all_samp_corrs <-  unlist(lapply(all_ds_corr_around, function(sub_data){
            lapply(sub_data, function(x) x[[paste0(corr_type)]])
          }))
        }
        stopifnot(!is.null(all_samp_corrs))
        all_samp_corrs <- na.omit(all_samp_corrs)  
        
        all_empPvals <- sapply(all_obs_corr, function(x) {
          (sum(all_samp_corrs >= x) + 1)/(length(all_samp_corrs) + 1)
        })
        
        stopifnot(file.path(hicds, exprds) %in% names(all_ds_corr_around))
        dsOnly <- all_ds_corr_around[[file.path(hicds, exprds)]]
        if(corr_type == "meanCorrLeftRight"){
          dsOnly_samp_corrs_left <- unlist(lapply(dsOnly, function(x) x[[paste0("meanCorr_left")]]))
          dsOnly_samp_corrs_right <- unlist(lapply(dsOnly, function(x) x[[paste0("meanCorr_right")]]))
          dsOnly_samp_corrs <- c(dsOnly_samp_corrs_left, dsOnly_samp_corrs_right)
        } else {
          dsOnly_samp_corrs <- unlist(lapply(dsOnly, function(x) x[[paste0(corr_type)]]))  
        }
        stopifnot(!is.null(dsOnly_samp_corrs))
        dsOnly_samp_corrs <- na.omit(dsOnly_samp_corrs)
    
        dsOnly_empPvals <- sapply(all_obs_corr, function(x) {
          (sum(dsOnly_samp_corrs >= x) + 1)/(length(dsOnly_samp_corrs) + 1)
        })
        dt <- data.frame(all_empPvals, dsOnly_empPvals)
        colnames(dt) <- c( paste0("empPval-", samp_type, "-", corr_type, " - allDS"),
                           paste0("empPval-", samp_type, "-", corr_type, " - DSonly"))
        
        dt
      } # end-iterating over type of corr
      corr_empPval_dt
    }
    stopifnot(rownames(corr_empPval_dt) == names(all_obs_corr))
    
    
    outFile <- file.path(outFold, hicds, exprds, "all_cors_empPval_dt.Rdata")
    dir.create(dirname(outFile), recursive = TRUE)
    save(all_cors_empPval_dt, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))

  } # end iterating over exprds
} # end iterating over hicds


#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
# all_ds_corr_around_TADs_sameKb_file <- file.path("COEXPR_AROUND_TADS", "sameKb", "all_ds_around_TADs_corr.Rdata")
# all_ds_corr_around_TADs_fixKb_file <- file.path("COEXPR_AROUND_TADS", "fixKb", fixKbSize, "all_ds_around_TADs_corr.Rdata")
# all_ds_corr_around_TADs_sameNbr_file <- file.path("COEXPR_AROUND_TADS", "sameNbr", "all_ds_around_TADs_corr.Rdata")
# stopifnot(file.exists(all_ds_corr_around_TADs_sameKb_file))
# stopifnot(file.exists(all_ds_corr_around_TADs_fixKb_file))
# stopifnot(file.exists(all_ds_corr_around_TADs_sameNbr_file))
# all_ds_corr_around_TADs_sameKb <- eval(parse(text = load(all_ds_corr_around_TADs_sameKb_file)))
# all_ds_corr_around_TADs_fixKb <- eval(parse(text = load(all_ds_corr_around_TADs_fixKb_file)))
# all_ds_corr_around_TADs_sameNbr <- eval(parse(text = load(all_ds_corr_around_TADs_sameNbr_file)))
# 
# ### sameKb - meanCorr
# all_meanCorr_sameKb <- unlist(lapply(all_ds_corr_around_TADs_sameKb, function(sub_data){
#   lapply(sub_data, function(x) x[["meanCorr"]])
# }))
# all_meanCorr_sameKb <- na.omit(all_meanCorr_sameKb)
# 
# ### sameKb - meanCorrRight
# all_meanCorrRight_sameKb <- unlist(lapply(all_ds_corr_around_TADs_sameKb, function(sub_data){
#   lapply(sub_data, function(x) x[["meanCorr_right"]])
# }))
# all_meanCorrRight_sameKb <- na.omit(all_meanCorrRight_sameKb)
# 
# ### sameKb - meanCorr
# all_meanCorrLeft_sameKb <- unlist(lapply(all_ds_corr_around_TADs_sameKb, function(sub_data){
#   lapply(sub_data, function(x) x[["meanCorr_left"]])
# }))
# all_meanCorrLeft_sameKb <- na.omit(all_meanCorrLeft_sameKb)
# 
# 
# all_empPvalsLeft_sameKb <- sapply(all_obs_corr, function(x) {
#   (sum(all_meanCorrLeft_sameKb >= x) + 1)/length(all_meanCorrLeft_sameKb)
# })
# 

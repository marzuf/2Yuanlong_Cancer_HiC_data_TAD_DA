
# Rscript check_script_pipeline.R

script_name <- "check_script_pipeline.R"

startTime <- Sys.time()

require(foreach)
require(doMC)
registerDoMC(40)


mainFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
# check if in pipeline region list

hicds = all_hicds[1]
exprds = all_exprds[[paste0(hicds)]]


script5_name <- "5sameNbr_runPermutationsCorr"
script5_v0folder <- "CREATE_SAMPLE_AROUND_TADS_SAMENBR_TWOSIDED"

script7_name <- "7sameNbr_runPermutationsMeanTADCorr"
script7_v0folder <- file.path("COEXPR_AROUND_TADS", "sameNbr")

script10_name <- "10sameNbr_runEmpPvalMeanTADCorr"
script10_v0folder <- "SAMPLE_MEANCORR_EMPPVALS"

script11_name <- "11sameNbr_runEmpPvalCombined"
script11_v0folder <- file.path("CREATE_EMPPVALCOMB", "sameNbr")

script19_name <- "19sameNbr_SAM_emp_measurement"
script19_v0folder <- file.path("SAM_EMP_MEASUREMENT_MEANCORR", "sameNbr")

# results for FC should be the same
script19_fc1_name <- "19onlyFC_SAM_emp_measurement"
script19_fc2_name <- "19onlyFCandCorr_SAM_emp_measurement"


foo <- foreach(hicds = all_hicds) %dopar% {
  foo <- foreach(exprds=all_exprds[[paste0(hicds)]]) %do% {
    
    cat("> START: ", hicds, " - ", exprds, "\n")
    
    ### CHECK SCRIPT 5 -> PERMUTATIONS
    cat("... check script 5\n")
    pip_file <- file.path(pipFolder, hicds, exprds, script5_name, "sample_around_TADs_sameNbr.Rdata")
    v0_file <- file.path(script5_v0folder, hicds, exprds, "sample_around_TADs_sameNbr.Rdata" )
    stopifnot(file.exists(pip_file))
    stopifnot(file.exists(v0_file))
    pip_data <- get(load(pip_file))
    v0_data <- get(load(v0_file))
    stopifnot(all.equal(pip_data, v0_data))
    
    ### CHECK SCRIPT 7 -> MEAN CORR
    cat("... check script 7\n")
    pip_file <- file.path(pipFolder, hicds, exprds, script7_name, "meanCorr_sample_around_TADs_sameNbr.Rdata")
    v0_file <- file.path(script7_v0folder, "all_ds_around_TADs_corr.Rdata" )
    stopifnot(file.exists(pip_file))
    stopifnot(file.exists(v0_file))
    pip_data <- get(load(pip_file))
    v0_data <- get(load(v0_file))
    pip_data2 <- unlist(lapply(pip_data, function(x) x[["meanCorr"]]))
    v0_data2 <- unlist(lapply(v0_data[[file.path(hicds, exprds)]], function(x)x[["meanCorr"]]))
    stopifnot(all.equal(pip_data2, v0_data2))
    
    ### CHECK SCRIPT 10 -> MEAN TAD CORR FOR THE PERMUT
    cat("... check script 10\n")
    pip_file <- file.path(pipFolder, hicds, exprds, script10_name, "emp_pval_meanCorr.Rdata")
    v0_file <- file.path(script10_v0folder, hicds, exprds, "all_cors_empPval_dt.Rdata" )
    stopifnot(file.exists(pip_file))
    stopifnot(file.exists(v0_file))
    pip_data <- get(load(pip_file))
    v0_data <- get(load(v0_file))
    v0_data2 <- setNames(v0_data[,"empPval-sameNbr-meanCorr - allDS"], rownames(v0_data))
    stopifnot(all.equal(pip_data, v0_data2))
    
    ### CHECK SCRIPT 11 -> COMBINED PVALS
    cat("... check script 11\n")
    pip_file <- file.path(pipFolder, hicds, exprds, script11_name, "emp_pval_combined.Rdata")
    v0_file <- file.path(script11_v0folder, hicds, exprds, "meanCorr_meanLogFC_notAdjCombEmpPval.Rdata" )
    stopifnot(file.exists(pip_file))
    stopifnot(file.exists(v0_file))
    pip_data <- get(load(pip_file))
    v0_data <- get(load(v0_file))
    stopifnot(all.equal(pip_data, v0_data))
    
    ### CHECK SCRIPT 19 -> SAM
    cat("... check script 19\n")
    pip_file <- file.path(pipFolder, hicds, exprds, script19_name, "meanCorr_empFDR.Rdata")
    v0_file <- file.path(script19_v0folder, hicds, exprds, "all_empFDR.Rdata" )
    stopifnot(file.exists(pip_file))
    stopifnot(file.exists(v0_file))
    pip_data <- get(load(pip_file))
    v0_data <- get(load(v0_file))
    v0_data2 <- unlist(v0_data["sample_meanCorr_allDS"], recursive = FALSE, use.names = FALSE)
    stopifnot(all.equal(pip_data[[1]], v0_data2[[1]]))
    stopifnot(all.equal(pip_data[[2]], v0_data2[[2]]))
    
    
    ### CHECK SCRIPT 19 FC -> SAM
    cat("... check script 19 FC\n")
    pip_file <- file.path(pipFolder, hicds, exprds, script19_fc1_name, "empFDR_list.Rdata")
    v0_file <- file.path(pipFolder, hicds, exprds, script19_fc2_name, "empFDR_list.Rdata")
    stopifnot(file.exists(pip_file))
    stopifnot(file.exists(v0_file))
    pip_data <- get(load(pip_file))
    v0_data <- get(load(v0_file))
    v0_data2 <- v0_data[grepl("logFC", names(v0_data))]
    all.equal(v0_data2, pip_data)
    
    
    
    
  }
}
  
  


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




  

  
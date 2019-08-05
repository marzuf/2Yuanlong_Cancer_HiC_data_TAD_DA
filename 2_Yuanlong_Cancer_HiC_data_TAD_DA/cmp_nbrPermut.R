options(scipen=100)

# Rscript cmp_nbrPermut.R

startTime <- Sys.time()

script_name <- "cmp_nbrPermut.R"

startTime <- Sys.time()

cat("> START ", script_name," \n")

setDir <- "/media/electron"
setDir <- ""

require(foreach)
require(doMC)

registerDoMC(40)

combPvalFolder_1 <- file.path("CREATE_EMPPVALCOMB", "sameNbr")
combPvalFolder_0 <- file.path("CREATE_EMPPVALCOMB_10000permut", "sameNbr")

script9_name_1 <- "9_runEmpPvalMeanTADLogFC"
script9_name_0 <- "910000_runEmpPvalMeanTADLogFC"

finalTable_0 <- file.path("CREATE_FINAL_TABLE_10000permut", "all_result_dt.Rdata")
finalTable_1 <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
finalTable_1_txt <- file.path("CREATE_FINAL_TABLE", "all_result_dt.txt")


CMP_SCRIPT9 <- FALSE
CMP_CREATETABLE <- TRUE
CMP_COMBPVAL <- TRUE

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

pipFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE","OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))

outFolder <- file.path("CMP_NBRPERMUT")
dir.create(outFolder, recursive=TRUE)


all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))


combPvalFolder <- file.path("CREATE_EMPPVALCOMB", "sameNbr")
stopifnot(dir.exists(combPvalFolder))

if(CMP_SCRIPT9) {
  hicds="GSE105381_HepG2_40kb"
  all_comp_permut_DT <- foreach(hicds = all_hicds, .combine='rbind') %do% {
    exprds="TCGAlihc_norm_lihc"  
    ds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat("... start ", hicds, " - ", exprds, "\n")
      
      file_v0 <- file.path(pipFolder, hicds, exprds, script9_name_0, "emp_pval_meanLogFC.Rdata")
      stopifnot(file.exists(file_v0))
      empPvalFC_0 <- eval(parse(text = load(file_v0)))
      
      file_v1 <- file.path(pipFolder, hicds, exprds, script9_name_1, "emp_pval_meanLogFC.Rdata")
      stopifnot(file.exists(file_v1))
      empPvalFC_1 <- eval(parse(text = load(file_v1)))
      
      stopifnot(setequal(names(empPvalFC_1), names(empPvalFC_0)))
      stopifnot(length(empPvalFC_0) == length(empPvalFC_1))
      
      data.frame(
        hicds=hicds,
        exprds=exprds,
        region=names(empPvalFC_0),
        empPvalFC_0 = empPvalFC_0,
        empPvalFC_1 = empPvalFC_1[names(empPvalFC_0)],
        stringsAsFactors = FALSE
      )
    } # end-foreach iterating exprds
    ds_dt
  } # end-foreach iterating over all_comp_permut_DT
  
  totDS <- length(unique(paste0(all_comp_permut_DT$hicds, "_",all_comp_permut_DT$exprds)))
  
  myx <- all_comp_permut_DT[,"empPvalFC_0"]
  myy <- all_comp_permut_DT[,"empPvalFC_1"]
  
  myxlab <- "empPvalFC 10^4 permut"
  myylab <- "empPvalFC 10^5 permut"
  
  outFile <- file.path(outFolder, paste0("empPvalFC_", "permut0", "_vs_", "permut1", ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x=myx,
    y=myy,
    xlab = myxlab,
    ylab = myylab,
    main = paste0("10^5 vs. 10^4 permut")
  )
  mtext(side=3, text = paste0("allDS - n=", totDS))
  addCorr(x = myx, y = myy, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
} # end-if comp script 9


if(CMP_COMBPVAL) {
  hicds="GSE105381_HepG2_40kb"
  all_comp_permut_DT <- foreach(hicds = all_hicds, .combine='rbind') %do% {
    exprds="TCGAlihc_norm_lihc"  
    ds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat("... start ", hicds, " - ", exprds, "\n")
      
      file_v0 <- file.path(combPvalFolder_0, hicds, exprds, "meanCorr_meanLogFC_notAdjCombEmpPval.Rdata")
      stopifnot(file.exists(file_v0))
      empPvalComb_0 <- eval(parse(text = load(file_v0)))
      
      file_v1 <- file.path(combPvalFolder_1, hicds, exprds, "meanCorr_meanLogFC_notAdjCombEmpPval.Rdata")
      stopifnot(file.exists(file_v1))
      empPvalComb_1 <- eval(parse(text = load(file_v1)))
      
      stopifnot(setequal(names(empPvalComb_1), names(empPvalComb_0)))
      stopifnot(length(empPvalComb_0) == length(empPvalComb_1))
      
      data.frame(
        hicds=hicds,
        exprds=exprds,
        region=names(empPvalComb_0),
        empPvalComb_0 = empPvalComb_0,
        empPvalComb_1 = empPvalComb_1[names(empPvalComb_0)],
        stringsAsFactors = FALSE
      )
    } # end-foreach iterating exprds
    ds_dt
  } # end-foreach iterating over all_comp_permut_DT
  
  totDS <- length(unique(paste0(all_comp_permut_DT$hicds, "_",all_comp_permut_DT$exprds)))
  
  myx <- all_comp_permut_DT[,"empPvalComb_0"]
  myy <- all_comp_permut_DT[,"empPvalComb_1"]
  
  myxlab <- "empPvalComb 10^4 permut"
  myylab <- "empPvalComb 10^5 permut"
  
  outFile <- file.path(outFolder, paste0("empPvalComb_", "permut0", "_vs_", "permut1", ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x=myx,
    y=myy,
    xlab = myxlab,
    ylab = myylab,
    main = paste0("10^5 vs. 10^4 permut")
  )
  mtext(side=3, text = paste0("allDS - n=", totDS))
  addCorr(x = myx, y = myy, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
} # end-if comp script 9


if(CMP_CREATETABLE) {
  hicds="GSE105381_HepG2_40kb"
  all_comp_permut_DT <- foreach(hicds = all_hicds, .combine='rbind') %do% {
    exprds="TCGAlihc_norm_lihc"  
    ds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat("... start ", hicds, " - ", exprds, "\n")
      
      stopifnot(file.exists(finalTable_0))
      final_DT_0 <- eval(parse(text = load(finalTable_0)))
      
      stopifnot(file.exists(finalTable_1))
      final_DT_1 <- eval(parse(text = load(finalTable_1)))
      
      # final_DT_1 <- read.delim(finalTable_1_txt, header=T, stringsAsFactors = FALSE)
      
      final_DT_0 <- final_DT_0[order(final_DT_0$hicds, final_DT_0$exprds, final_DT_0$region),]
      final_DT_1 <- final_DT_1[order(final_DT_1$hicds, final_DT_1$exprds, final_DT_1$region),]
      
      # colnames(final_DT_0) <- paste0(colnames(final_DT_0), "_0")
      # colnames(final_DT_1) <- paste0(colnames(final_DT_1), "_1")
      
      
      all_vars <- c("meanLogFC_thresh_FDR0.1", "meanLogFC_thresh_FDR0.2",
                    "meanCorr_thresh_FDR0.1", "meanCorr_thresh_FDR0.2")
      
      same_vars <- c("hicds", "exprds", "region")
      #"region_genes", "meanLogFC", 
       #              "meanCorr", "ratioDown")
      
      for(curr_var in same_vars) {
        stopifnot(final_DT_1[,curr_var] == final_DT_0[,curr_var])
      }
      all_vars %in% colnames(final_DT_0)
      all_vars %in% colnames(final_DT_1)
      final_DT <- merge(final_DT_0[,c(same_vars, all_vars)], 
                        final_DT_1[,c(same_vars, all_vars)], 
                        by =same_vars, suffixes = c("_0","_1"))
     
      
      final_DT
      
    } # end-foreach iterating exprds
    ds_dt
  } # end-foreach iterating over all_comp_permut_DT
  
  save(all_comp_permut_DT, "all_comp_permut_DT_ft.Rdata")
  
  totDS <- length(unique(paste0(all_comp_permut_DT$hicds, "_",all_comp_permut_DT$exprds)))
  
  
  for(vars in all_vars) {
    
    myx <- all_comp_permut_DT[,paste0(vars, "_0")]
    myy <- all_comp_permut_DT[,paste0(vars, "_1")]
    
    myxlab <- paste0(vars, " 10^4 permut")
    myylab <- paste0(vars, " 10^5 permut")
    
    outFile <- file.path(outFolder, paste0(vars, "permut0", "_vs_", "permut1", ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x=myx,
      y=myy,
      xlab = myxlab,
      ylab = myylab,
      main = paste0("10^5 vs. 10^4 permut")
    )
    mtext(side=3, text = paste0("allDS - n=", totDS))
    addCorr(x = myx, y = myy, bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
  }
  
} # end-if comp create table



##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




options(scipen=100)

# Rscript boundary_effect_TAD.R   # to run all datasets in one shot

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "boundary_effect_TAD.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

tiesMeth <- "min"

buildData <- TRUE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

script0_name <- "0_prepGeneData"
script4_name <- "4_runMeanTADCorr"


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4


pipFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

corrPvalsFolder <- file.path("SAMPLE_MEANCORR_EMPPVALS")
stopifnot(dir.exists(corrPvalsFolder))

outFolder <- "BOUNDARY_EFFECT_TAD"
dir.create(outFolder, recursive=TRUE)

corr_type <- "meanCorr"
samp_type <- "sameNbr"
permut_type <- "allDS"
fixSizeKb <- ""
sampName <- paste0("empPval-", samp_type, "-", corr_type, " - ", permut_type)



pvalFolder <- "CREATE_EMPPVALCOMB"
samp_type <- "sameNbr"
fixSizeKb <- ""
stopifnot(dir.exists(pvalFolder))



samp_file <- file.path("COEXPR_AROUND_TADS", samp_type, fixSizeKb, "all_ds_around_TADs_corr.Rdata")
stopifnot(file.exists(samp_file))
all_ds_corr_around <- eval(parse(text = load(samp_file)))

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


if(buildData) {
  
  ### BUILD SIGNIF ALONG FDR THRESH
  cat("... start building signif. along FDR thresh\n")
  all_corr_tad_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds_corr_tad_DT <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      obs_corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
      stopifnot(file.exists(obs_corr_file))
      all_obs_corr <- eval(parse(text = load(obs_corr_file)))
      
      ds_corr_around <- all_ds_corr_around[[file.path(hicds, exprds)]]
      
      all_samp_corrs <-  unlist(lapply(ds_corr_around, function(x){
               x[[paste0("meanCorr")]]
      }))        
                
      
      all_regs <- names(all_obs_corr)
      stopifnot(all_regs %in% names(all_samp_corrs))
      

      
      combPvalFile <- file.path(pvalFolder, samp_type, fixSizeKb, hicds, exprds, "meanCorr_meanLogFC_notAdjCombEmpPval.Rdata")
      stopifnot(file.exists(combPvalFile))
      tad_combPvals <- eval(parse(text=load(combPvalFile)))
      adj_tad_combPvals <- p.adjust(tad_combPvals, method="BH")
      
      stopifnot(all_regs %in% names(adj_tad_combPvals))
      
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        region = all_regs,
        obs_meanCorr = all_obs_corr[all_regs],
        around_meanCorr = all_samp_corrs[all_regs],
        tad_adjPvalComb = adj_tad_combPvals[all_regs],
        stringsAsFactors = FALSE
      )          
    } # end iterating over exprds
  } # end iterating over hicds
  
  outFile <- file.path(outFolder, "all_corr_tad_dt.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(all_corr_tad_dt, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else { # end-if build data
  outFile <- file.path(outFolder, "all_corr_tad_dt.Rdata")  
  all_corr_tad_dt <- eval(parse(text = load(outFile)))
}

all_corr_tad_dt <- all_corr_tad_dt[all_corr_tad_dt$obs_meanCorr > 0 & all_corr_tad_dt$around_meanCorr > 0,]

all_corr_tad_dt$corrRatio <- all_corr_tad_dt$obs_meanCorr/all_corr_tad_dt$around_meanCorr

all_corr_tad_dt$corrRatio_log2 <- log2(all_corr_tad_dt$obs_meanCorr/all_corr_tad_dt$around_meanCorr)

all_corr_tad_dt$tad_adjPvalComb_log10 <- -log10(all_corr_tad_dt$tad_adjPvalComb)
      
      
all_corr_tad_dt$subtype_col <- all_cols[paste0(all_cmps[paste0(all_corr_tad_dt$exprds)])]
all_corr_tad_dt$dataset <- paste0( all_corr_tad_dt$hicds, " - ", all_corr_tad_dt$exprds)
all_datasets <- unique(all_corr_tad_dt$dataset) 
nDS <- length(unique( all_datasets ))
      
###################################
################################### PLOT EACH DATASET SEPARATELY
###################################
x_var = "tad_adjPvalComb_log10"
y_var = "corrRatio_log2"

foo <- foreach(ds = all_datasets) %dopar% {
  
  sub_dt <- all_corr_tad_dt[all_corr_tad_dt$dataset == ds,]
  
  myx <- sub_dt[,paste0(x_var)]
  myy <- sub_dt[,paste0(y_var)]
  outFile <- file.path(outFolder, paste0(gsub(" - ", "_", ds), "_", y_var, "_vs_", x_var, "_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x = myx,
    y = myy,
    cex.axis = axisCex,
    cex.lab = axisCex,
    xlab = paste0(gsub("_", " ", gsub("_log10", "", x_var), " [-log10]")),
    ylab =  paste0(gsub("_", " " , gsub("_log2", "", y_var)), " [>0, log2]"),
    main = paste0(y_var, " vs. ", x_var)
  )
  mtext(side=3, text = paste0(ds, " (n=", nrow(sub_dt),")"), font=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}


###################################
################################### PLOT FOR ALL DATASETS
###################################

      
x_var = "tad_adjPvalComb_log10"
y_var = "corrRatio_log2"
myx <- all_corr_tad_dt[,paste0(x_var)]
myy <- all_corr_tad_dt[,paste0(y_var)]
outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = myx,
  y = myy,
  cex.axis = axisCex,
  cex.lab = axisCex,
  xlab = paste0(gsub("_", " ", gsub("_log10", "", x_var)) , " [-log10]"),
  ylab = paste0(gsub("_", " " , gsub("_log2", "", y_var)), " [>0, log2]"),
  main = paste0(y_var, " vs. ", x_var)
)
mtext(side=3, text = paste0(nDS, " DS (n=", nrow(all_corr_tad_dt),")"), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


x_var = "tad_adjPvalComb_log10"
y_var = "corrRatio_log2"
myx <- all_corr_tad_dt[,paste0(x_var)]
myy <- all_corr_tad_dt[,paste0(y_var)]
outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_plotColType.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x = myx,
  y = myy,
  cex.axis = axisCex,
  cex.lab = axisCex,
  xlab = paste0(gsub("_", " ", gsub("_log10", "", x_var)), " [-log10]"),
  ylab = paste0(gsub("_", " " , gsub("_log2", "", y_var)), " [>0, log2]"),
  main = paste0(y_var, " vs. ", x_var),
  col = all_corr_tad_dt$subtype_col,
  pch=16,
  cex=0.7
)
addSubtypeLeg(mypos="bottomright", bty="n")
mtext(side=3, text = paste0(nDS, " DS (n=", nrow(all_corr_tad_dt),")"), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

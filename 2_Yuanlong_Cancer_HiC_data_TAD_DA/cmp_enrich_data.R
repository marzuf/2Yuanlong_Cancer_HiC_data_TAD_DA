options(scipen=100)

# Rscript cmp_enrich_data.R

startTime <- Sys.time()

script_name <- "cmp_enrich_data.R"

startTime <- Sys.time()

cat("> START ", script_name," \n")

outFolder <- "CMP_ENRICH_DATA"
dir.create(outFolder, recursive = TRUE)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4


all_hmm <- c("K562", "Hepg2", "Gm12878")
i=1;j=2
for(i in 1:(length(all_hmm) -1)) {
  
  i_hmm <- all_hmm[i]
  i_file <- file.path("ENHANCER_ENRICH", i_hmm, "all_chromHMM_tadCount_DT.Rdata")
  stopifnot(file.exists(i_file))
  i_DT <- eval(parse(text = load(i_file)))
  
  all_vars <- colnames(i_DT)
  
  totDS <- length(unique(paste0(i_DT$hicds, "_", i_DT$exprds)))
  
  
  for(j in (i+1):length(all_hmm) ) {
    
    j_hmm <- all_hmm[j]
    j_file <- file.path("ENHANCER_ENRICH", j_hmm, "all_chromHMM_tadCount_DT.Rdata")
    stopifnot(file.exists(j_file))
    j_DT <- eval(parse(text = load(j_file)))
    stopifnot(setequal(colnames(i_DT), colnames(j_DT)))
    
    same_vars <- c("region", "hicds", "exprds", "region_adjPvalComb")
    all_vars <- all_vars[!all_vars %in% same_vars]
    
    for(vars in same_vars) {
      stopifnot(i_DT[,vars] == j_DT[,vars])
    }
    
    for(vars in all_vars) {
      
      i_data <- i_DT[,vars]
      j_data <- j_DT[,vars]
      
      
      outFile <- file.path(outFolder, paste0("allDS", "_", gsub("\\(|\\)|\\/","_", gsub(" ", "", vars)), "_", i_hmm, "_vs_", j_hmm, "densplot.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      densplot(
        x=i_data,
        y=j_data,
        xlab = i_hmm,
        ylab = j_hmm,
        main = paste0(gsub("\\(|\\)|\\/","_", gsub(" ", "", vars)))
      )
      mtext(side=3, text = paste0("allDS - n=", totDS, " - ", i_hmm, " vs. ", j_hmm))
      addCorr(x = i_data, y = j_data, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      
      
      
    }
    
    
    
  }
  
}




##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




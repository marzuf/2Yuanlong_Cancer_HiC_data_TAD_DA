options(scipen=100)

# Rscript combined_pvals.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript combined_pvals.R   # to run all datasets in one shot

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "combined_pvals.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
require(ggplot2)
require(reshape2)


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("utils_fct.R")

script0_name <- "0_prepGeneData"
script9_name <- "9_runEmpPvalMeanTADLogFC"
# ../Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/9_runEmpPvalMeanTADLogFC/emp_pval_meanLogFC.Rdata

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4

myHeightGG <- 7
myWidthGG <- myHeightGG 


pval_thresh <- 0.05


pipOutFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))


corrPvalsFolder <- file.path("SAMPLE_MEANCORR_EMPPVALS")
stopifnot(dir.exists(corrPvalsFolder))

outFolder <- "COMBINED_PVALS"
dir.create(outFolder, recursive=TRUE)


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




for(hicds in all_hicds) {
  
  for(exprds in all_exprds[[paste0(hicds)]]) {
    # SAMPLE_MEANCORR_EMPPVALS/Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/all_cors_empPval_dt.Rdata
    
    
    empPval_logFC_file <- file.path(pipOutFolder, hicds, exprds, script9_name, "emp_pval_meanLogFC.Rdata" )
    stopifnot(file.exists(empPval_logFC_file))
    empPval_logFC <- eval(parse(text = load(paste0(empPval_logFC_file))))
    
    all_regs <- names(empPval_logFC)
    
    
    empPval_meanCorr_file <- file.path(corrPvalsFolder, hicds, exprds, "all_cors_empPval_dt.Rdata")
    stopifnot(file.exists(empPval_meanCorr_file))
    empPval_meanCorr_dt <- eval(parse(text = load(paste0(empPval_meanCorr_file))))
    
    stopifnot(setequal(rownames(empPval_meanCorr_dt), all_regs))
    
    empPval_meanCorr_dt <- empPval_meanCorr_dt[all_regs,]
    
    adj_combinedPvals_dt <- apply(empPval_meanCorr_dt, 2, function(corrType_col) {
      stopifnot(length(corrType_col) == length(empPval_logFC))
      combined_pvals <- sapply(1:length(corrType_col), function(i) {
        pval_logFC <- empPval_logFC[i]
        stopifnot(is.numeric(pval_logFC))
        pval_meanCorr <- corrType_col[i]
        stopifnot(is.numeric(pval_meanCorr))
        stouffer(c(pval_logFC, pval_meanCorr), two.tails=TRUE)
      })
      combined_pvals <- p.adjust(combined_pvals, method = "BH")   #### ==> !!! the saved pvals are adjusted !!!
      names(combined_pvals) <- names(corrType_col)
      combined_pvals
    })
    # JUST CHECK WITH THE 1st COLUMN THAT I DID WHAT I WANTED    
    check_meanCorrPvals <- setNames(empPval_meanCorr_dt[,1], rownames(empPval_meanCorr_dt))
    stopifnot(names(check_meanCorrPvals) == names(empPval_logFC))
    check_pvalsComb <- sapply(1:length(empPval_logFC), function(x) {
      stouffer(c(empPval_logFC[x], check_meanCorrPvals[x]), two.tails=TRUE)
    })
    check_pvalsComb <- p.adjust(check_pvalsComb, method="BH")
    stopifnot(check_pvalsComb == adj_combinedPvals_dt[,1])
    
    colnames(adj_combinedPvals_dt) <- gsub("empPval", "adjComb_empPval", colnames(adj_combinedPvals_dt))
    all_adjCombinedPvals_dt <- adj_combinedPvals_dt
    outFile <- file.path(outFolder, hicds, exprds, "all_adjCombinedPvals_dt.Rdata")
    dir.create(dirname(outFile), recursive = TRUE)
    save(all_adjCombinedPvals_dt, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    ###################################################################################
    ################################################################################### DENSITY CORR PVALS
    ###################################################################################
    corrPvals_list1 <- apply(empPval_meanCorr_dt, 2, as.list)
    corrPvals_list <- lapply(corrPvals_list1, as.numeric)
    corrPvals_list_log10 <- lapply(corrPvals_list, function(x) -log10(as.numeric(x)))
    
    outFile <- file.path(outFolder, hicds, exprds, paste0("allCorr_corrEmpPval_densplot.", plotType))
    dir.create(dirname(outFile), recursive = TRUE)
    do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
    plot_multiDens(
      corrPvals_list_log10,
      plotTit = paste0(hicds, " - ", exprds),
      my_xlab = paste0("emp. p-val [-log10]")
                   )
    mtext(side=3, text=paste0("all types"))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    toKeep1 <- which(grepl("meanCorr - allDS", names(corrPvals_list_log10)))
    toKeep2 <- which(grepl("meanCorrLeftRight - allDS", names(corrPvals_list_log10)))
    
    outFile <- file.path(outFolder, hicds, exprds, paste0("meanCorr_meanCorrLeftRight_allDS_corrEmpPval_densplot.", plotType))
    dir.create(dirname(outFile), recursive = TRUE)
    do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
    plot_multiDens(
      corrPvals_list_log10[c(toKeep1, toKeep2)],
      plotTit = paste0(hicds, " - ", exprds),
      my_xlab = paste0("emp. p-val [-log10]")
    )
    mtext(side=3, text=paste0("allDS - meanCorr and meanCorrLeftRight"))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    ###################################################################################
    ################################################################################### DENSITY COMBINED PVALS
    ###################################################################################
    
    combinedPvals_list1 <- apply(adj_combinedPvals_dt, 2, as.list)
    combinedPvals_list <- lapply(combinedPvals_list1, as.numeric)
    combinedPvals_list_log10 <- lapply(combinedPvals_list1, function(x) -log10(as.numeric(x)))
    
    outFile <- file.path(outFolder, hicds, exprds, paste0("allCorr_adjCombEmpPval_densplot.", plotType))
    dir.create(dirname(outFile), recursive = TRUE)
    do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
    plot_multiDens(
      combinedPvals_list_log10,
      plotTit = paste0(hicds, " - ", exprds),
      my_xlab = paste0("emp. p-val [-log10]")
    )
    mtext(side=3, text=paste0("all types"))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    toKeep1 <- which(grepl("meanCorr - allDS", names(combinedPvals_list_log10)))
    toKeep2 <- which(grepl("meanCorrLeftRight - allDS", names(combinedPvals_list_log10)))
    
    outFile <- file.path(outFolder, hicds, exprds, paste0("meanCorr_meanCorrLeftRight_allDS_adjCombEmpPval_densplot.", plotType))
    dir.create(dirname(outFile), recursive = TRUE)
    do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
    plot_multiDens(
      combinedPvals_list_log10[c(toKeep1, toKeep2)],
      plotTit = paste0(hicds, " - ", exprds),
      my_xlab = paste0("emp. p-val [-log10]")
    )
    mtext(side=3, text=paste0("allDS - meanCorr and meanCorrLeftRight"))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    ###################################################################################
    ################################################################################### BOXPLOT NBR SIGNIF
    ###################################################################################
    
    melt_comb_dt <- melt(all_adjCombinedPvals_dt)
    colnames(melt_comb_dt) <- c("region", "pvalType", "pval")
    stopifnot(is.numeric(melt_comb_dt$pval))
    stopifnot(setequal(melt_comb_dt$region, all_regs))
    
    
    nSignif <- c(by(melt_comb_dt, melt_comb_dt$pvalType, function(sub_dt) sum(sub_dt$pval <= pval_thresh)))
    
    nSignif_dt <- data.frame(
      corrType = names(nSignif),
      nSignif = as.numeric(nSignif),
      stringsAsFactors = FALSE
    )
    
    nSignif_dt$samp_type <- gsub("adjComb_empPval-(.+)-mean.+", "\\1", nSignif_dt$corrType)
    nSignif_dt$corr_type <- gsub("adjComb_empPval-.+-(mean.+)", "\\1", nSignif_dt$corrType)
    nSignif_dt$corr_type <- gsub(" - ", "\n", nSignif_dt$corr_type)
    
    p_nsignif <-  ggplot(nSignif_dt, aes(x = samp_type, y = nSignif, fill = samp_type)) + 
      geom_bar(position="dodge", stat="identity") +
      facet_grid(~corr_type, switch="x") + 
      coord_cartesian(expand = FALSE) +
      ggtitle(paste0(hicds, " - ", exprds), subtitle = paste0("# signif. TADs - adj. comb. p-val. thresh = ", pval_thresh) )+
      scale_x_discrete(name="")+
      scale_y_continuous(name=paste0("# signif. TADs"),
                         breaks = scales::pretty_breaks(n = 20))+
      # scale_fill_brewer(palette="YlOrRd")+
      labs(fill  = "") +
      theme( # Increase size of axis lines
        strip.text = element_text(size = 10),
        # top, right, bottom and left
        # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", size=16),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.minor.y = element_line(colour = "grey"),
        strip.text.x = element_text(size = 10),
        axis.line.x = element_line(size = .2, color = "black"),
        axis.line.y = element_line(size = .3, color = "black"),
        axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        legend.background =  element_rect(),
        legend.key = element_blank(),
        legend.title = element_text(face="bold")
      )
    outFile <- file.path(outFolder, hicds, exprds, paste0("nSignif_by_corrType_sampType_signifThresh", pval_thresh , "_barplot.", plotType))
    dir.create(dirname(outFile), recursive = TRUE)
    ggsave(p_nsignif, filename = outFile, height = myHeightGG, width = myWidthGG*1.4)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
  } # end - for iterating over exprds
} # end - for iterating over hicds









##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



    
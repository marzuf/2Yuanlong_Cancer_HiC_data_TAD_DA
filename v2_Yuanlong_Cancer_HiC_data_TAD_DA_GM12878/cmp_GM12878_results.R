# ^

require(foreach)
require(doMC)
registerDoMC(40)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- file.path("CMP_GM12878_RESULTS")
dir.create(outFolder)

corMet <- "pearson"

plotType <- "png"
myHeight <- myWidth <- 400

main_dt <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA//GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
gm_dt <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))

main_dt$ds <- file.path(main_dt$hicds, main_dt$exprds)
gm_dt$ds <- file.path(gm_dt$hicds, gm_dt$exprds)


main_dt$tad_adjCombPval_log10 <- -log10(main_dt$tad_adjCombPval)
gm_dt$tad_adjCombPval_log10 <- -log10(gm_dt$tad_adjCombPval)

all_dt <- merge(main_dt, gm_dt, by=c("exprds", "entrezID"), suffixes = c("_MAIN", "_GM"))


gm12878_yl_dt <- merge(main_dt, gm_dt, by=c("exprds", "entrezID"), suffixes = c("_yl", "_gm12878"))
save(gm12878_yl_dt, file=file.path(outFolder, "gm12878_yl_dt.Rdata"), version=2)


all_vars <- c("gene_rank", "tad_rank", "tad_adjCombPval", "adj.P.Val", "tad_adjCombPval_log10")
curr_var = all_vars[1]

all_corr_dt <- foreach(curr_var = all_vars, .combine='rbind')%dopar% {
  
  my_x <- all_dt[,paste0(curr_var, "_MAIN")]
  my_y <- all_dt[,paste0(curr_var, "_GM")]
  
  outFile <- file.path(outFolder, paste0(curr_var, "_GM12878_vs_MAIN_densplot.png"))
  do.call("png", list(outFile, height=myHeight, width=myWidth))
    densplot(x = my_x, y=my_y,
             main=paste0(curr_var),
             xlab="MAIN", ylab="GM12878") 
    mtext(side=3, text=paste0("# DS = ", length(unique(file.path(all_dt$hicds_GM,all_dt$hicds_MAIN, all_dt$exprds))), 
                              " (# MAIN = ", length(unique(file.path(main_dt$ds))), 
                              "; # GM12878 = ",length(unique(file.path(gm_dt$ds))), ")" ))
    addCorr(x=my_x, y=my_y, legPos = "topleft", bty="n")
    foo <- dev.off()
    
    var_corr_dt <- foreach(ds = unique(main_dt$ds), .combine='rbind') %dopar% {
      
      sub_main_dt <- main_dt[main_dt$ds == ds,]
      sub_all_dt <- merge(sub_main_dt, gm_dt, by=c("exprds", "entrezID"), suffixes = c("_MAIN", "_GM")) 
      cor_val <- cor(sub_all_dt[,paste0(curr_var, "_MAIN")],
                     sub_all_dt[,paste0(curr_var, "_GM")], method=corMet
      )
      data.frame(
        dataset = ds,
        variable = curr_var, 
        corrCoeff = cor_val,
        stringsAsFactors = FALSE
      )
    }
    
    var_corr_dt <- var_corr_dt[order(var_corr_dt$corrCoeff, decreasing = TRUE),]
    outFile <- file.path(outFolder, paste0(curr_var, "_GM12878_vs_MAIN_corrCoeff_barplot.svg"))
    do.call("svg", list(outFile, height=5, width=5))
    barplot(var_corr_dt$corrCoeff,
            xlab=paste0("all datasets (n=", nrow(var_corr_dt), ")"), ylab=paste0(corMet, "'s corr. coeff."),
            main=paste0(curr_var, " - corr. MAIN and GM12878"))
    foo <- dev.off()
    
    var_corr_dt
    
}
save(all_corr_dt, file=file.path(outFolder, "all_corr_dt.Rdata"), version=2)
# load("CMP_GM12878_RESULTS/all_corr_dt.Rdata)

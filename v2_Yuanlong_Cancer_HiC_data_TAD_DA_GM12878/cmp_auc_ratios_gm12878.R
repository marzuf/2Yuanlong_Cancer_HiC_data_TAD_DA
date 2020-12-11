dt1 <- get(load("../MANUSCRIPT_FIGURES/FIG_1/FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata"))
dt2 <- get(load("FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata"))

all_dt <- merge(dt1, dt2, by="exprds", suffixes = c("_hicds", "_gm12878"), all=T)

outFolder <- "CMP_AUC_RATIOS_GM12878"
dir.create(outFolder)

auc_gm12878_yl_dt <- merge(dt1, dt2, by="exprds", suffixes = c("_yl", "_gm12878"), all=T)
save(auc_gm12878_yl_dt, file=file.path(outFolder, "auc_gm12878_yl_dt.Rdata"), version=2)


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")



xlab <- "all hicds"
ylab <- "GM12878"

myx <- all_dt$fcc_auc_hicds
myy <- all_dt$fcc_auc_gm12878

plotType <- "png"
outFile <- file.path(outFolder, paste0("cmp_fcc_auc.", plotType))
do.call(plotType, list(outFile, height=400, width=400))
plot(myx,myy, xlab=xlab, ylab=ylab, pch=16, cex=0.7, main="YL data")
addCorr(x=myx,y=myy, legPos="topleft", bty="n")
mtext(side=3, text=paste0("nDS = ", nrow(all_dt)))
curve(1*x, add=T,col="grey")

dev.off()


# Rscript ctcf_curves_positiveControl.R

library("readxl")
library(doMC)
library(foreach)
library(stringr)


require(ggpubr)
require(ggsci)
registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
plotType <- "png"
myHeight <- myWidth <- 400
plotCex <- 1.2

plotTypeGG <- "svg"
ggHeight <- 6
ggWidth <- 5

fontFamily <- "Hershey"


runFolder <- "." 

tad_dt <- get(load("CONVERT_BD_TADS/all_tads_dt.Rdata"))  

source("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/ctcf_da_utils.R")

nBreaks <- 100
step_breaks <- 1/nBreaks


buildTable <- F

outFolder <- file.path("CTCF_CURVES_POSITIVECONTROL")
dir.create(outFolder, recursive = TRUE)

init_ctcf_dt <- read_excel("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/13059_2020_2108_MOESM2_ESM.xlsx", sheet="CTCFs")
init_ctcf_dt <- as.data.frame(init_ctcf_dt)
init_ctcf_dt <- init_ctcf_dt[, 1:7]
init_ctcf_dt$chr <- as.character(init_ctcf_dt$chr)

if(buildTable){
  # assign ctcf BS to tads
  ctcf_dt <- init_ctcf_dt[init_ctcf_dt$chr %in% tad_dt$chromo,]
  stopifnot(nrow(ctcf_dt) > 0)
  stopifnot(is.numeric(ctcf_dt$start))
  stopifnot(is.numeric(ctcf_dt$end))
  stopifnot(ctcf_dt$start <= ctcf_dt$end)
  ctcf_dt$region <- NA
  
  i=1
  ctcf2tad_dt <- foreach(i = 1:nrow(ctcf_dt), .combine='rbind') %dopar% {
    
    chr <- ctcf_dt$chr[i]
    ctcf_start <- ctcf_dt$start[i]
    ctcf_end <- ctcf_dt$end[i]
    
    subtad_dt <- tad_dt[tad_dt$chromo == chr,]
    stopifnot(nrow(subtad_dt) > 0)
    
    # assign if start after tad start and end before tad end
    test1 <- which(ctcf_start >= subtad_dt$start & ctcf_end <= subtad_dt$end)
    test2 <- which(subtad_dt$start <= ctcf_start & subtad_dt$end >= ctcf_end)
    stopifnot(test1==test2)
    stopifnot(length(test1) == 0 | length(test1) == 1)
    if(length(test1) == 1) {
      ctcf_dt$region[i] <- subtad_dt$region[test1]
    } else {
      ctcf_dt$region[i] <- NA
    }
    ctcf_dt[i,]  
  }
  
  outFile <- file.path(outFolder, "ctcf2tad_dt.Rdata")
  save(ctcf2tad_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile,"\n"))
  
}else {
  # outFile="CTCF_CURVES_POSITIVECONTROL/ctcf2tad_dt.Rdata"
  outFile <- file.path(outFolder, "ctcf2tad_dt.Rdata")
  ctcf2tad_dt <- get(load(outFile))
}

# stop("-ok\n")
# ctcf2tad_dt=get(load("CTCF_AND_DA_ALLDS/ctcf2tad_dt.Rdata"))

colnames(tad_dt)[colnames(tad_dt) == "start"] <- "tad_start"
colnames(tad_dt)[colnames(tad_dt) == "end"] <- "tad_end"

both_dt <- merge(tad_dt, ctcf2tad_dt, by =c("region"), all=FALSE)

stopifnot(nrow(both_dt) > 0)
both_dt$ctcf_midpos <- (both_dt$start + both_dt$end)/2
stopifnot(both_dt$end <= both_dt$tad_end)
stopifnot(both_dt$start >= both_dt$tad_start)

both_dt$relative_position <- (both_dt$ctcf_midpos - both_dt$tad_start)/(both_dt$tad_end - both_dt$tad_start)
stopifnot(both_dt$relative_position >= 0)
stopifnot(both_dt$relative_position <= 1)

# ctcf2tad_dt$relative_position <- (ctcf2tad_dt$ctcf_midpos - both_dt$tad_start)/(both_dt$tad_end - both_dt$tad_start)
# stopifnot(both_dt$relative_position >= 0)
# stopifnot(both_dt$relative_position <= 1)
# 
# # rel_pos_breaks <- seq(from=0, to = 1, length.out=nBreaks)
# # rel_pos_labs <- get_fract_lab0(vect_values=both_dt$relative_position, range_levels = rel_pos_breaks)
# rel_pos_levels <- gsub("^<=0$", "0",get_level_labs(rel_pos_breaks))
# both_dt$rel_pos_lab <- rel_pos_labs

both_dt$rel_pos_lab <- both_dt$relative_position %/% step_breaks
both_dt$rel_pos_lab <- factor(both_dt$rel_pos_lab, levels=0:nBreaks)
stopifnot(!is.na(both_dt$rel_pos_lab))

my_theme <-theme_bw()+
  theme(
    axis.text.x = element_blank(),
    legend.position="top",
    # panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    # panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    # panel.background = element_rect(fill = "transparent"),
    # panel.grid.major.x =  element_blank(),
    # panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    plot.title = element_text(hjust=0.5, size = 14, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.text = element_text(size=12, hjust=0.5, vjust=0.5)
  )

### corrected in v3

curr_col <- "relative_position"
curr_fun1 <- "length"
curr_fun2 <- "sum"

agg_byTAD_dt <- aggregate(as.formula(paste0(curr_col, " ~ region + rel_pos_lab")), 
                          data = both_dt, FUN=curr_fun1)

col_leg <- ""
# paste0(names(table(agg_byTAD_dt$signif_lab)), "\n", as.numeric(table(agg_byTAD_dt$signif_lab)))

# myTit <- paste0("# DS = ", length(unique(file.path(agg_byTAD_dt$hicds,agg_byTAD_dt$exprds))), 
#                 "; # TADs = ", length(unique(file.path(agg_byTAD_dt$hicds,agg_byTAD_dt$exprds, agg_byTAD_dt$region ))))
myTit <- paste0("# TADs = ", length(unique(file.path(agg_byTAD_dt$region ))))

agg_all_dt <- aggregate(as.formula(paste0(curr_col, " ~ rel_pos_lab")), data=agg_byTAD_dt, FUN = curr_fun2)

# tmp_agg_dt <- unique(agg_byTAD_dt[,c("hicds", "exprds", "region", "signif_lab")])
# col_leg <- paste0(names(table(tmp_agg_dt$signif_lab)), "\n", as.numeric(table(tmp_agg_dt$signif_lab)))

agg_all_dt[,paste0(curr_col, "_rescaled")] <-  agg_all_dt[,paste0(curr_col)]/length(unique(agg_byTAD_dt$region))

p1 <-  ggplot(data = agg_all_dt, aes_string(x="rel_pos_lab", y = paste0(curr_col, ""), group=1))+#, color="signif_lab", group="signif_lab")) + 
  labs(x="relative position in TAD", y=paste0(curr_fun1, " agg. ", curr_col))+
  # scale_color_d3(labels = col_leg)+
  ggtitle(myTit) +
  geom_line() + 
  geom_smooth(method = "loess")+
  scale_y_continuous(breaks=scales::pretty_breaks(n = 10)) +
  my_theme
  
  outFile <- file.path(outFolder, paste0(curr_col, "_aggBy_", curr_fun1, "_along_relPosInTAD_bySignif_lineplot.", plotTypeGG))
  ggsave(p1, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
  
  p1 <- ggplot(data = agg_all_dt, aes_string(x="rel_pos_lab", y = paste0(curr_col, "_rescaled"), group=1))+#, color="signif_lab", group="signif_lab")) + 
    labs(x="relative position in TAD", y=paste0(curr_fun1, " agg. ", curr_col, " (resc.)"))+
    # scale_color_d3(labels = col_leg)+
    ggtitle(myTit) +
    geom_line() + 
    geom_smooth(method = "loess")+
    scale_y_continuous(breaks=scales::pretty_breaks(n = 10)) +
    my_theme
  

outFile <- file.path(outFolder, paste0(curr_col, "_aggBy_", curr_fun1, "_along_relPosInTAD_bySignif_lineplot_rescaled.", plotTypeGG))
ggsave(p1, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))


##################################################################################### by orientation
########### by orientation

agg_byTAD_dt <- aggregate(as.formula(paste0(curr_col, " ~ region + rel_pos_lab + orientation")), data = both_dt, FUN=curr_fun1)

col_leg <- ""

myTit <- paste0("# TADs = ", length(unique(file.path(agg_byTAD_dt$region ))))

agg_all_dt <- aggregate(as.formula(paste0(curr_col, " ~ rel_pos_lab + orientation")), data=agg_byTAD_dt, FUN = curr_fun2)

tmp_agg_dt <- unique(agg_byTAD_dt[,c( "region", "orientation")])
col_leg <- paste0(names(table(tmp_agg_dt$orientation)), "\n", as.numeric(table(tmp_agg_dt$orientation)))

agg_all_dt[,paste0(curr_col, "_rescaled")] <-  agg_all_dt[,paste0(curr_col)]/table(tmp_agg_dt$orientation)[agg_all_dt$orientation] 


p2 <- ggplot(data = agg_all_dt, aes_string(x="rel_pos_lab", y = paste0(curr_col), color="orientation", group="orientation")) + 
  labs(x="relative position in TAD", y=paste0(curr_fun1, " agg. ", curr_col), color="")+
  ggtitle(myTit) +
  scale_color_d3(labels = col_leg)+
  geom_line() + 
  geom_smooth(method = "loess")+
  scale_y_continuous(breaks=scales::pretty_breaks(n = 10))+
  my_theme
outFile <- file.path(outFolder, paste0(curr_col, "_aggBy_", curr_fun1, "_along_relPosInTAD_byOrientation_lineplot.", plotTypeGG))
ggsave(p2, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))

p2 <- ggplot(data = agg_all_dt, aes_string(x="rel_pos_lab", y = paste0(curr_col, "_rescaled"), color="orientation", group="orientation")) + 
  labs(x="relative position in TAD", y=paste0(curr_fun1, " agg. ", curr_col, " (resc.)"), color="")+
  ggtitle(myTit) +
  scale_color_d3(labels = col_leg)+
  geom_line() + 
  geom_smooth(method = "loess")+
  scale_y_continuous(breaks=scales::pretty_breaks(n = 10))+
  my_theme
outFile <- file.path(outFolder, paste0(curr_col, "_aggBy_", curr_fun1, "_along_relPosInTAD_byOrientation_lineplot_rescaled.", plotTypeGG))
ggsave(p2, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))





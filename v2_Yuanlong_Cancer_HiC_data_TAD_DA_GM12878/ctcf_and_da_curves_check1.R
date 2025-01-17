
# Rscript ctcf_and_da_curves_check1.R


################# AJOUTER LE # TOT DANS LA LEGENDE !!!

source("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/ctcf_da_utils.R")

library(ggplot2)
library(ggsci)

pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
all_obs_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_obs_exprds <- sapply(all_obs_hicds, function(x) list.files(file.path(pipFolder, x)))


all_obs_exprds <- "TCGAluad_mutKRAS_mutEGFR"
all_obs_hicds <- "GM12878_40kb"
# inFolder <- "CTCF_AND_DA_ALLDS_CHECK1"
inFolder <- "CTCF_AND_DA_ALLDS"


nBreaks <- 100
step_breaks <- 1/nBreaks

plotTypeGG <- "svg"
ggHeight <- 5
ggWidth <- 6

outFolder <- file.path("CTCF_AND_DA_CURVES_CHECK1")
dir.create(outFolder, recursive = TRUE)

inFile <- file.path(inFolder, "ctcf2tad_dt.Rdata")
ctcf2tad_dt <- get(load(inFile))
# stopifnot(ctcf2tad_dt$hicds %in% all_obs_hicds)

inFile <- file.path("PSEUDO_FINAL_TABLE/all_result_dt.Rdata")
final_dt <- get(load(inFile))
ds_final_dt <- final_dt
ds_final_dt <- final_dt[final_dt$hicds %in% unlist(all_obs_hicds) & final_dt$exprds %in% unlist(all_obs_exprds),  ]
colnames(ds_final_dt)[colnames(ds_final_dt) == "start"] <- "tad_start"
colnames(ds_final_dt)[colnames(ds_final_dt) == "end"] <- "tad_end"


both_dt <- merge(ds_final_dt, ctcf2tad_dt, by =c("hicds", "region"), all=FALSE)
stopifnot(nrow(both_dt) > 0)
both_dt$ctcf_midpos <- (both_dt$start + both_dt$end)/2
stopifnot(both_dt$end <= both_dt$tad_end)
stopifnot(both_dt$start >= both_dt$tad_start)

both_dt$relative_position <- (both_dt$ctcf_midpos - both_dt$tad_start)/(both_dt$tad_end - both_dt$tad_start)
stopifnot(both_dt$relative_position >= 0)
stopifnot(both_dt$relative_position <= 1)

# rel_pos_breaks <- seq(from=0, to = 1, length.out=nBreaks)
# rel_pos_labs <- get_fract_lab0(vect_values=both_dt$relative_position, range_levels = rel_pos_breaks)
# rel_pos_levels <- gsub("^<=0$", "0",get_level_labs(rel_pos_breaks))
# both_dt$rel_pos_lab <- rel_pos_labs

both_dt$rel_pos_lab <- both_dt$relative_position %/% step_breaks
both_dt$rel_pos_lab <- factor(both_dt$rel_pos_lab, levels=0:nBreaks)
stopifnot(!is.na(both_dt$rel_pos_lab))


head(both_dt)


both_dt$rel_pos_lab2 <- both_dt$relative_position


############################################################################ by signif

all_to_plot <- list(
  c("relative_position", "length"), 
  c("MotifScore", "mean"), 
  c("MotifScore", "sum"), 
  c("ChipSeqScore", "mean"), 
  c("ChipSeqScore", "sum")
)

pthresh <- 0.01

both_dt$signif_lab <- ifelse(both_dt$adjPvalComb <= pthresh, "signif.", "not signif.")

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

  
for(i in 1:length(all_to_plot)) {
  
  
  curr_col <- all_to_plot[[i]][[1]]
  curr_fun <- all_to_plot[[i]][[2]]
  
  ########### by signif
  
  agg_byTAD_dt <- aggregate(as.formula(paste0(curr_col, " ~ hicds + exprds + region + rel_pos_lab + signif_lab")), data = both_dt, FUN=curr_fun)
  
  
  # ggplot(both_dt, aes(x = "relative_position", col = "signif_lab")) + geom_density()
  # 
  # ggdensity(both_dt,
  #           x = "relative_position",
  #           y = "..density..",
  #           # combine = TRUE,                  # Combine the 3 plots
  #           xlab = "relative_position",
  #           # add = "median",                  # Add median line.
  #           rug = FALSE,                      # Add marginal rug
  #           color = "signif_lab",
  #           fill = "signif_lab",
  #           palette = "jco"
  # ) 
  # both_dt2 <- both_dt
  # both_dt2$rel_pos_lab <- as.numeric(both_dt2$rel_pos_lab)
  # ggdensity(both_dt2,
  #           x = "rel_pos_lab",
  #           y = "..density..",
  #           # combine = TRUE,                  # Combine the 3 plots
  #           xlab = "relative_position",
  #           # add = "median",                  # Add median line.
  #           rug = FALSE,                      # Add marginal rug
  #           color = "signif_lab",
  #           fill = "signif_lab",
  #           palette = "jco"
  # ) 
    
  
  myTit <- paste0("# DS = ", length(unique(file.path(agg_byTAD_dt$hicds,agg_byTAD_dt$exprds))), 
                  "; # TADs = ", length(unique(file.path(agg_byTAD_dt$hicds,agg_byTAD_dt$exprds, agg_byTAD_dt$region ))))
  
  agg_all_dt <- aggregate(as.formula(paste0(curr_col, " ~ rel_pos_lab + signif_lab")), data=agg_byTAD_dt, FUN = "mean")
  
  tmp_agg_dt <- unique(agg_byTAD_dt[,c("hicds", "exprds", "region", "signif_lab")])
  col_leg <- paste0(names(table(tmp_agg_dt$signif_lab)), "\n", as.numeric(table(tmp_agg_dt$signif_lab)))
  
  
  p1 <- ggplot(data = agg_all_dt, aes_string(x="rel_pos_lab", y = paste0(curr_col), color="signif_lab", group="signif_lab")) + 
    labs(x="relative position in TAD", y=paste0(curr_fun, " agg. ", curr_col), color="")+
    scale_color_d3(labels = col_leg)+
    ggtitle(myTit) +
    geom_line() + 
    geom_smooth(method = "loess")+
    scale_y_continuous(breaks=scales::pretty_breaks(n = 10)) +
    my_theme
  
  outFile <- file.path(outFolder, paste0(curr_col, "_aggBy_", curr_fun, "_along_relPosInTAD_bySignif_lineplot.", plotTypeGG))
  ggsave(p1, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
  
  # agg_dt2 <- aggregate(as.formula(paste0("relative_position", " ~ rel_pos_lab + signif_lab")), data = both_dt, FUN="length")
  # ggplot(data = agg_dt2, aes_string(x="rel_pos_lab", y = paste0("relative_position"), color="signif_lab", group="signif_lab")) + 
  #   labs(x="relative position in TAD", y=paste0(curr_fun, " agg. ", curr_col), color="")+
  #   scale_color_d3(labels = col_leg)+
  #   ggtitle(myTit) +
  #   geom_line() + 
  #   geom_smooth(method = "loess")+
  #   scale_y_continuous(breaks=scales::pretty_breaks(n = 10)) +
  #   my_theme
  
  
  
  
  
  ########### by triplet class
  
  agg_byTAD_dt <- aggregate(as.formula(paste0(curr_col, " ~ hicds + exprds + region + rel_pos_lab + Triplet_class")), data = both_dt, FUN=curr_fun)
  
  col_leg <- paste0(names(table(agg_byTAD_dt$Triplet_class)), "\n", as.numeric(table(agg_byTAD_dt$Triplet_class)))
  

  myTit <- paste0("# DS = ", length(unique(file.path(agg_byTAD_dt$hicds,agg_byTAD_dt$exprds))), 
                  "; # TADs = ", length(unique(file.path(agg_byTAD_dt$hicds,agg_byTAD_dt$exprds, agg_byTAD_dt$region ))))
  
  agg_all_dt <- aggregate(as.formula(paste0(curr_col, " ~ rel_pos_lab + Triplet_class")), data=agg_byTAD_dt, FUN = "mean")
  
  p2 <- ggplot(data = agg_all_dt, aes_string(x="rel_pos_lab", y = paste0(curr_col), color="Triplet_class", group="Triplet_class")) + 
    labs(x="relative position in TAD", y=paste0(curr_fun, " agg. ", curr_col), color="")+
    ggtitle(myTit) +
    scale_color_d3(labels = col_leg)+
    geom_line() + 
    geom_smooth(method = "loess")+
    scale_y_continuous(breaks=scales::pretty_breaks(n = 10))+
    my_theme
  
  
  outFile <- file.path(outFolder, paste0(curr_col, "_aggBy_", curr_fun, "_along_relPosInTAD_byTripletClass_lineplot.", plotTypeGG))
  ggsave(p2, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
  ########### by orientation
  
  agg_byTAD_dt <- aggregate(as.formula(paste0(curr_col, " ~ hicds + exprds + region + rel_pos_lab + orientation")), data = both_dt, FUN=curr_fun)
  
  col_leg <- paste0(names(table(agg_byTAD_dt$orientation)), "\n", as.numeric(table(agg_byTAD_dt$orientation)))
  
  
  myTit <- paste0("# DS = ", length(unique(file.path(agg_byTAD_dt$hicds,agg_byTAD_dt$exprds))), 
                  "; # TADs = ", length(unique(file.path(agg_byTAD_dt$hicds,agg_byTAD_dt$exprds, agg_byTAD_dt$region ))))
  
  agg_all_dt <- aggregate(as.formula(paste0(curr_col, " ~ rel_pos_lab + orientation")), data=agg_byTAD_dt, FUN = "mean")
  
  p2 <- ggplot(data = agg_all_dt, aes_string(x="rel_pos_lab", y = paste0(curr_col), color="orientation", group="orientation")) + 
    labs(x="relative position in TAD", y=paste0(curr_fun, " agg. ", curr_col), color="")+
    ggtitle(myTit) +
    scale_color_d3(labels = col_leg)+
    geom_line() + 
    geom_smooth(method = "loess")+
    scale_y_continuous(breaks=scales::pretty_breaks(n = 10))+
    my_theme
  
  
  outFile <- file.path(outFolder, paste0(curr_col, "_aggBy_", curr_fun, "_along_relPosInTAD_byOrientation_lineplot.", plotTypeGG))
  ggsave(p2, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
  
}




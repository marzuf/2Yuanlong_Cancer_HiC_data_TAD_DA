

# Rscript boxplot_FDRthresh_nbrSignif.R

options(scipen = 100)

script_name <- "boxplot_FDRthresh_nbrSignif.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()


require(reshape2)
require(ggplot2)
require(hrbrthemes)
require(ggthemes)
require(foreach)
require(doMC)
registerDoMC(40)

pipOutFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "OUTPUT_FOLDER")
inFolder <- file.path("MEANCORR_MEANFC_FDR")

outFolder <- "BOXPLOT_FDRTHRESH_NBRSIGNIF"
dir.create(outFolder, recursive = TRUE)

plotType <- "png"
myHeightGG  <- 7
myWidthGG <- 12

args <- commandArgs(trailingOnly = TRUE)
hicds <- args[1]
exprds <- args[2]


if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}


hicds = all_hicds[1]
all_dt <- foreach(hicds = all_hicds, .combine='rbind') %do% {
  
  
  
  exprds = all_exprds[[paste0(hicds)]][1]
  all_dt_hicds <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    
    cat("... start building DT for :", hicds, " - ", exprds, "\n")
    
    curr_inFolder <- file.path(inFolder, hicds, exprds)
    
    v0_file <- file.path(curr_inFolder, "signifTADs_FC_and_CorrV0.Rdata")
    stopifnot(file.exists(v0_file))
    
    v0_FDR <- eval(parse(text = load(v0_file)))
    
    v0_DT <- do.call(rbind,lapply(v0_FDR, data.frame))
    v0_DT$FDR_cut_off <- rownames(v0_DT)
    rownames(v0_DT) <- NULL
    # colnames(v0_DT) <- paste0(colnames(v0_DT), "_v0")
    colnames(v0_DT) <- paste0(colnames(v0_DT), "_meanCorrV0")
    
    
    allcorr_file <- file.path(curr_inFolder, "signifTADs_FC_and_allCorrs.Rdata")
    stopifnot(file.exists(allcorr_file))
    
    allCorrs_FDR <- eval(parse(text = load(allcorr_file)))
    
    
    dt_list_allCorrs <- lapply(allCorrs_FDR, function(sublist1) lapply(sublist1, function(sublist2) do.call(rbind, lapply(sublist2,  data.frame))))
    stopifnot(length(dt_list_allCorrs) > 0)
    
    i=1
    allCorrs_DT <- foreach(i = 1:length(dt_list_allCorrs), .combine='cbind') %do% {
      dt_list <- dt_list_allCorrs[[i]]
      samp_type <- names(dt_list_allCorrs)[i]
      stopifnot(length(dt_list) > 0)
      j=1
      sub_dt <- foreach(j = 1:length(dt_list), .combine='cbind') %dopar% {
        tmp_dt <- dt_list[[j]]
        tmp_dt$FDR_cut_off <- rownames(tmp_dt)
        rownames(tmp_dt) <- NULL
        corr_type <- names(dt_list)[j]
        colnames(tmp_dt) <- paste0(colnames(tmp_dt), "_", samp_type, "_", corr_type)
        tmp_dt
      } 
      sub_dt
    }
    rownames(allCorrs_DT) <- NULL
    
    
    allCorrs_DT <- cbind(allCorrs_DT, v0_DT)
    
    stopifnot(length(unique(as.list(allCorrs_DT[,grepl("FDR_cut_off", colnames(allCorrs_DT))]))) == 1)
    allCorrs_DT$FDR_cut_off <- allCorrs_DT[,grep("FDR_cut_off", colnames(allCorrs_DT))[1]]
    stopifnot(length(unique(as.list(allCorrs_DT[,grepl("FDR_cut_off", colnames(allCorrs_DT))]))) == 1)
    
    stopifnot(length(unique(as.list(allCorrs_DT[,grepl("fc_cut_off", colnames(allCorrs_DT))]))) == 1)
    allCorrs_DT$fc_cut_off <- allCorrs_DT[,grep("fc_cut_off", colnames(allCorrs_DT))[1]]
    stopifnot(length(unique(as.list(allCorrs_DT[,grepl("fc_cut_off", colnames(allCorrs_DT))]))) == 1)
    
    
    allCorrs_DT <- allCorrs_DT[,!grepl("FDR_cut_off.+", colnames(allCorrs_DT))]
    allCorrs_DT <- allCorrs_DT[,!grepl("fc_cut_off.+", colnames(allCorrs_DT))]
    
    allCorrs_DT$hicds <- hicds
    allCorrs_DT$exprds <- exprds
    allCorrs_DT
  } # end-foreach iterating over exprds
  
  all_dt_hicds
  
} # end-foreach iterating over hicds

outFile <- file.path(outFolder, paste0("all_dt.Rdata"))
save(all_dt, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

plot_all_dt <- melt(all_dt, id = c("hicds", "exprds", "FDR_cut_off"))
plot_all_dt$variable <- as.character(plot_all_dt$variable)
plot_all_dt$value <- as.numeric(as.character(plot_all_dt$value))
stopifnot(!is.na(plot_all_dt$value))

all_vars <- c("nTADs_signif", "fc_cut_off", "corr_cut_off")

all_vars_tit <- setNames(c("# signif. TADs", "logFC cutoff", "meanCorr cutoff"), all_vars)

var_to_plot = all_vars[1]
var_to_plot = all_vars[2]

nDS <- length(unique(paste0(all_dt$hicds, "_", all_dt$exprds)))

for(var_to_plot in all_vars) {
  
  
  plot_dt <- plot_all_dt[grepl(paste0(var_to_plot), plot_all_dt$variable),]
  plot_dt$variable <- gsub(paste0(var_to_plot, "_"), "", plot_dt$variable) # will need to redo because nTADS written !!!
  plot_dt$variable <- gsub("_sample", "", plot_dt$variable)
  plot_dt$variable <- gsub("_meanCorr", "\nmeanCorr", plot_dt$variable)
  
  stopifnot("variable" %in% colnames(plot_dt))
  stopifnot("value" %in% colnames(plot_dt))
  
  p_var <-  ggplot(plot_dt, aes(x = variable, y = value, fill = variable)) + 
    geom_boxplot()+
    facet_grid(~FDR_cut_off, switch="x") + 
    coord_cartesian(expand = FALSE) +
    ggtitle(paste0(all_vars_tit[paste0(var_to_plot)], " by variable FDR cutoff"), subtitle = paste0(samp_type, " - nDS = ", nDS))+
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(all_vars_tit[paste0(var_to_plot)]),
                       breaks = scales::pretty_breaks(n = 20))+
    # scale_fill_brewer(palette="YlOrRd")+
    labs(fill  = "meanCorr type") +
    theme( # Increase size of axis lines
      strip.text = element_text(size = 12),
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
  
  outFile <- file.path(outFolder, paste0(var_to_plot, "_allDS_allSampTypes_boxplot.", plotType))
  ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  all_sampTypes <- unique(unlist(lapply(strsplit(x =plot_dt$variable, split="\n"), function(x)x[[1]])))
  samp_type=all_sampTypes[1]
  for(samp_type in all_sampTypes) {
    
    
    st_plot_dt <- plot_dt[grepl(samp_type, plot_dt$variable),]
    
    st_plot_dt$variable <- gsub(paste0(samp_type, "\n"), "", st_plot_dt$variable)
    
    nDS_sub <- length(unique(paste0(st_plot_dt$hicds, "_", st_plot_dt$exprds)))
    
    p_var_sampType <-  ggplot(st_plot_dt, aes(x = variable, y = value, fill = variable)) + 
      geom_boxplot()+
      facet_grid(~FDR_cut_off, switch="x") + 
      coord_cartesian(expand = FALSE) +
      ggtitle(paste0(all_vars_tit[paste0(var_to_plot)], " by variable FDR cutoff"), subtitle = paste0(samp_type, " - nDS = ", nDS_sub))+
      scale_x_discrete(name="")+
      scale_y_continuous(name=paste0(all_vars_tit[paste0(var_to_plot)]),
                         breaks = scales::pretty_breaks(n = 20))+
      # scale_fill_brewer(palette="YlOrRd")+
      labs(fill  = "meanCorr type") +
      theme( # Increase size of axis lines
        strip.text = element_text(size = 12),
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
    
    # x1 <- st_plot_dt[st_plot_dt$FDR_cut_off == 0.5 & st_plot_dt$variable == "meanCorrLeft_allDS",]
    # x2 <- st_plot_dt[st_plot_dt$FDR_cut_off == 0.5 & st_plot_dt$variable == "meanCorrLeftRight_allDS",]
    # > all(x1$hicds == x2$hicds)
    # [1] TRUE
    # > all(x1$exprds == x2$exprds)
    # [1] TRUE
    # > all(x1$value == x2$value)
    # [1] FALSE
    
    outFile <- file.path(outFolder, paste0(var_to_plot, "_allDS_", samp_type, "_boxplot.", plotType))
    ggsave(plot = p_var_sampType, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
  } # end-iterating over sampTypes
  
} # end-iterating over variables to plot

#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))



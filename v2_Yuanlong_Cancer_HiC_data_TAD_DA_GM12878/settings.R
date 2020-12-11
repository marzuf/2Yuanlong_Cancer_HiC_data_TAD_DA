
nCpu <- 40

# V2
# labsymbol <- "\u25CF" # circle
# V3
labsymbol <- "\u25A0" # squares

source("full_dataset_names.R")

runFolder <-  "."

pipFolder <- file.path(runFolder, "PIPELINE/OUTPUT_FOLDER")

all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
all_obs_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_obs_exprds <- lapply(all_obs_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_obs_exprds) <- all_obs_hicds

get_exprds <- function(myfold) sapply(myfold, function(x) list.files(file.path(pipFolder, x)))

plotCex <- 1.4
axisCex <- 1.4
labCex <- 1.4
mainCex <- 1.4
subCex <- 1.2

myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 400, 7)

myHeightGG <- 7
myWidthGG <- 7

step0_folder <- "0_prepGeneData"
step8fcc_folder <- "8cOnlyFCC_runAllDown"

fontFamily <- "Hershey"

# ggsci package palettes:
#pal_npg()(100) # 10
#pal_aaas()(100) # 10
#pal_nejm()(100) # 8
#pal_lancet()(100) # 9
#pal_jama()(100) # 7
#pal_jco()(100) # 10
#pal_ucscgb()(100) # 26
#pal_d3()(100) # 10
#pal_locuszoom()(100) # 7
#pal_igv()(100)  # 51
#pal_cosmic()(100) # NF
#pal_uchicago()(100) # 9
#pal_startrek()(100) # 7
#pal_tron()(100) # 7
#pal_futurama()(100) # 12
#pal_rickandmorty()(100) # 12
#pal_simpsons()(100) # 16
#pal_gsea()(100) # 12
#pal_material()(100) # 10
require(ggsci)
observ_col <- pal_lancet()(3)[1]
permut_col <- pal_lancet()(3)[2]

tad_signif_col <- "dodgerblue3"
gene_signif_col <- "firebrick3"

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
# V2:
#all_cols[all_cols == "red"] <- "brown3"  # wt vs mut
#all_cols[all_cols == "blue"] <- "darkblue" # norm vs tum
#all_cols[all_cols == "green"] <- "forestgreen" # subtypes
# V3:
#all_cols[all_cols == "red"] <- "brown3"  # wt vs mut
#all_cols[all_cols == "blue"] <- "darkblue" # norm vs tum
#all_cols[all_cols == "green"] <- "forestgreen" # subtypes


#all_cols[all_cols == "red"] <- "violetred" #"chocolate"  # wt vs mut
#all_cols[all_cols == "blue"] <- "slateblue" # norm vs tum
#all_cols[all_cols == "green"] <- "slategray" # "yellow3" # subtypes

#all_cols[all_cols == "red"] <- "firebrick3" #"chocolate"  # wt vs mut
#all_cols[all_cols == "blue"] <- "navy" # norm vs tum
#all_cols[all_cols == "green"] <- "gray20" # "yellow3" # subtypes
all_cols[all_cols == "red"] <- "firebrick3" #"chocolate"  # wt vs mut
all_cols[all_cols == "blue"] <- "navy" # norm vs tum
all_cols[all_cols == "green"] <- "gray50" # "yellow3" # subtypes


options(save.defaults = list(version=2), scipen=100)

geneSignifThresh <- 0.01
tadSignifThresh <- 0.01

# FOR FIGURE 1:
require(ggplot2)
my_box_theme <- theme(
	text = element_text(family=fontFamily),
  panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.background = element_rect(fill = "transparent"),
  panel.grid.major.x =  element_blank(),
  panel.grid.minor.x =  element_blank(),
  axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
  axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
  plot.title = element_text(hjust=0.5, size = 16, face="bold"),
  plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
  legend.title = element_text(face="bold")
) 
  
require(scales)
noZero_breaks <- function (n = 5, ...) {
  scales:::force_all(n, ...)
  function(x) {
    breaks <- pretty(x, n, ...)
    breaks <- breaks[breaks > 0]
    names(breaks) <- c(attr(breaks, "labels"))
    c(1,breaks)
  }
}


exprBox_cond1Col <- "black"#"gray8"
exprBox_cond2Col <- "gray48"




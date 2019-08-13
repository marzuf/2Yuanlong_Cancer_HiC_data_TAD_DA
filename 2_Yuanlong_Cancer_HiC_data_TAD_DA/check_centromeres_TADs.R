options(scipen=100)


# Rscript check_centromeres_TADs.R

script_name <- "check_centromeres_TADs.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(rCGH)
require(foreach)
require(doMC)
require(lattice)

registerDoMC(40)

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

centroStart <- setNames(hg19$centromerStart, paste0("chr", hg19$chrom))
stopifnot(is.numeric(centroStart))
centroEnd <- setNames(hg19$centromerEnd, paste0("chr", hg19$chrom))
stopifnot(is.numeric(centroEnd))

mainFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

outFolder <- "CHECK_CENTROMERES_TADS"
dir.create(outFolder, recursive = TRUE)

hicds = all_hicds[1]
all_ctr_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt")
  stopifnot(file.exists(hicds_file))
  
  tadpos_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
  tadpos_dt <- tadpos_dt[grepl("_TAD", tadpos_dt$region),]
  stopifnot(nrow(tadpos_dt) > 0 )
  
  all_chr <- unique(tadpos_dt$chromo)
  stopifnot(all_chr %in% names(centroStart))
  stopifnot(all_chr %in% names(centroEnd))
  
  ctr_chr <- sapply(all_chr, function(chr) {
    
    ctr_s <- as.numeric(centroStart[chr])
    ctr_e <- as.numeric(centroEnd[chr])
    stopifnot(!is.na(ctr_s))
    stopifnot(!is.na(ctr_e))
    ctr_mid <- (ctr_s + ctr_e)/2
    
    ctr_mid_TAD <- tadpos_dt$region[which(tadpos_dt$chromo == chr &
                                            tadpos_dt$start <= ctr_mid &
                                              tadpos_dt$end >= ctr_mid)]
    
    ctr_start_TAD <- tadpos_dt$region[which(tadpos_dt$chromo == chr &
                                            tadpos_dt$start <= ctr_s &
                                            tadpos_dt$end >= ctr_s)]
    
    ctr_end_TAD <- tadpos_dt$region[which(tadpos_dt$chromo == chr &
                                            tadpos_dt$start <= ctr_e &
                                            tadpos_dt$end >= ctr_e)]
    list(ctr_mid_TAD=ctr_mid_TAD, ctr_start_TAD=ctr_start_TAD, ctr_end_TAD=ctr_end_TAD)
    
  })
  
  out_dt <- data.frame(t(ctr_chr))
  out_dt$chr <- rownames(out_dt)
  out_dt$hicds <- hicds
  out_dt <- out_dt[,c("hicds", "chr", "ctr_mid_TAD", "ctr_start_TAD", "ctr_end_TAD")]
  out_dt
}
outFile <- file.path(outFolder, "all_ctr_dt.Rdata")
save(all_ctr_dt, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

#######################################################################################################################################
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
# check if in pipeline region list

hicds = all_hicds[1]
exprds = all_exprds[[paste0(hicds)]]

all_ctrInRegion_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds_ctr_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    script0_name <- "0_prepGeneData"
    regionList_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
    stopifnot(file.exists(regionList_file))
    regionList <- eval(parse(text = load(regionList_file)))
    
    ctr_dt <- all_ctr_dt[all_ctr_dt$hicds == hicds,]
    
    ctr_dt$mid_TAD_inRegionList <- ctr_dt$ctr_mid_TAD %in% regionList
    ctr_dt$start_TAD_inRegionList <- ctr_dt$ctr_start_TAD %in% regionList
    ctr_dt$end_TAD_inRegionList <- ctr_dt$ctr_end_TAD %in% regionList
    
    all_cols <- colnames(ctr_dt)[colnames(ctr_dt) != "hicds"]
    
    ctr_dt$exprds <- exprds
    
    ctr_dt <- ctr_dt[,c("hicds", "exprds", all_cols)]
    ctr_dt
  }
  exprds_ctr_dt
}

outFile <- file.path(outFolder, "all_ctrInRegion_dt.Rdata")
save(all_ctrInRegion_dt, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))




#######################################################################################################################################
all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt")
  stopifnot(file.exists(hicds_file))
  
  tadpos_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
  tadpos_dt <- tadpos_dt[grepl("_TAD", tadpos_dt$region),]
  stopifnot(nrow(tadpos_dt) > 0 )
  
  all_cols <- colnames(tadpos_dt)
  tadpos_dt$hicds <- hicds
  tadpos_dt <- tadpos_dt[, c("hicds", all_cols)]
  
  tadpos_dt
  
}
outFile <- file.path(outFolder, "all_dt.Rdata")
save(all_dt, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
####################################################################################################################################### TAD size

all_dt$size <- all_dt$end - all_dt$start + 1

all_dt$size_log10 <- log10(all_dt$size)

#### PLOT THE SIGNIF PVAL FEATURES => compare dist Pval signif. vs. not signif.

# length(unique(all_dt$hicds)) # 18

nRows <- 5
nCols <- 4

all_vars <- c("size", "size_log10")
plot_var=all_vars[1]
for(plot_var in all_vars) {
  outFile <- file.path(outFolder, paste0("all_hicds_dist_TAD_", plot_var, "_density_lattice.", plotType))
  do.call(plotType, list(outFile, height=myHeight*3, width=myWidth*3))
  # myplot <- densityplot( formula(paste0("~",plot_var)), groups = hicds, data = all_dt, #auto.key = TRUE,
  myplot <- densityplot( formula(paste0("~",plot_var, "| hicds")), data = all_dt, #auto.key = TRUE,
                         # par.strip.text=list(cex=1), # width of the strip bar
                         par.strip.text = list(cex = 0.8, font = 4, col = "brown"),
                         layout = c(nCols, nRows), # column,row
                         scales=list(y=list(relation="free"),
                                     x=list(relation="free")
                         ),
                         main = paste0(plot_var))
  print(myplot)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}


####################################################################################################################################### # TADs by chromo

all_dt_chr <- aggregate(region ~ hicds + chromo, FUN=length, data = all_dt)

colnames(all_dt_chr)[colnames(all_dt_chr) == "region"] <- "chromo_nbrTADs"

all_vars <- c("chromo_nbrTADs")
plot_var=all_vars[1]
for(plot_var in all_vars) {
  outFile <- file.path(outFolder, paste0("all_hicds_dist_TAD_", plot_var, "_density_lattice.", plotType))
  do.call(plotType, list(outFile, height=myHeight*3, width=myWidth*3))
  # myplot <- densityplot( formula(paste0("~",plot_var)), groups = hicds, data = all_dt, #auto.key = TRUE,
  myplot <- densityplot( formula(paste0("~",plot_var, "| hicds")), data = all_dt_chr, #auto.key = TRUE,
                         # par.strip.text=list(cex=1), # width of the strip bar
                         par.strip.text = list(cex = 0.8, font = 4, col = "brown"),
                         layout = c(nCols, nRows), # column,row
                         scales=list(y=list(relation="free"),
                                     x=list(relation="free")
                         ),
                         main = paste0(plot_var))
  print(myplot)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}




#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

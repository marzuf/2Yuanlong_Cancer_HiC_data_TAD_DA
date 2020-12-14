require(doMC)
require(foreach)
registerDoMC(40)

# Rscript convert_bd_tads.R

outFolder <- "CONVERT_BD_TADS"
dir.create(outFolder, recursive = TRUE)

boundaries_dt <- read.delim("nanni_conserved_boundaries.csv", header=TRUE, sep=",",stringsAsFactors = FALSE)
conserved_dt <- boundaries_dt[boundaries_dt$Conservation >= 2,]
nrow(conserved_dt)
all_chromo <- unique(conserved_dt$chr)
chr = all_chromo[1]

all_tads_dt <- foreach(chr = all_chromo, .combine="rbind") %dopar% {
  
  chr_bd_dt <- conserved_dt[conserved_dt$chr == chr,]
  stopifnot(nrow(chr_bd_dt) > 0)
  
  max_size <- max(chr_bd_dt$end)
  
  # toxexample chr_bd_dt <- data.frame(start = c(10,60,100,120), end = c(40,80,120,150))
  # start end
  # 1    10  40
  # 2    60  80
  # 3   100 120
  # 4   120 150
  
  tad_dt <- data.frame(
    chromo = chr,
    start = c(0,chr_bd_dt$end[1:(nrow(chr_bd_dt)-1)]),
    end = c(chr_bd_dt$start), stringsAsFactors=FALSE)
  stopifnot(tad_dt$end >= tad_dt$start)
  tad_dt <- tad_dt[tad_dt$end != tad_dt$start,]
  tad_dt$region <- paste0(chr, "_TAD", 1:nrow(tad_dt))
  # start tad_end
  # 1     0      10
  # 2    40      60
  # 3    80     100
  stopifnot(max_size == (sum(chr_bd_dt$end-chr_bd_dt$start)+sum(tad_dt$end-tad_dt$start)))
  tad_dt
}
outFile <- file.path(outFolder, "all_tads_dt.Rdata")
save(all_tads_dt, version=2, file=outFile)
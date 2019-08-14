# Rscript cmp_filterCounts.R

script_name <- "cmp_filterCounts.R"

startTime <- Sys.time()

require(foreach)
require(doMC)
registerDoMC(40)


source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")

mainFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
# check if in pipeline region list

hicds = all_hicds[1]
exprds = all_exprds[[paste0(hicds)]]


script0_name <- "0_prepGeneData"
script0_name_v2 <- "0_prepGeneData_TEMP"

outFold <- "CMP_FILTERCOUNTS"
dir.create(outFold, recursive = TRUE)
pipLogFile <- file.path(outFold, "cmp_script0_minCounts_log_file.txt")
file.remove(pipLogFile)

foo <- foreach(hicds = all_hicds) %do% {
  foo <- foreach(exprds=all_exprds[[paste0(hicds)]]) %do% {
    
    txt <- paste0("> ", hicds, " - ", exprds, "\n")
    printAndLog(txt, pipLogFile)
    
    s0_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
    s0_v2_file <- file.path(pipFolder, hicds, exprds, script0_name_v2, "pipeline_geneList.Rdata")
    s0_genes <- get(load(s0_file))
    s2_genes <- get(load(s0_v2_file))
    txt <- paste0("# pipeline genes v0\t=\t", length(s0_genes), "\n")
    printAndLog(txt, pipLogFile)
    txt <- paste0("# pipeline genes v2\t=\t", length(s2_genes), "\n")
    printAndLog(txt, pipLogFile)
    txt <- paste0("# intersect pipeline genes\t=\t", length(intersect(s0_genes, s2_genes)), "\n")
    printAndLog(txt, pipLogFile)
    
    s0_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
    s0_v2_file <- file.path(pipFolder, hicds, exprds, script0_name_v2, "pipeline_regionList.Rdata")
    s0_regions <- get(load(s0_file))
    s2_regions <- get(load(s0_v2_file))
    txt <- paste0("# pipeline regions v0\t=\t", length(s0_regions), "\n")
    printAndLog(txt, pipLogFile)
    txt <- paste0("# pipeline regions v2\t=\t", length(s2_regions), "\n")
    printAndLog(txt, pipLogFile)
    txt <- paste0("# intersect pipeline regions\t=\t", length(intersect(s0_regions, s2_regions)), "\n")
    printAndLog(txt, pipLogFile)
    
  }
}

cat(paste0("... written: ", pipLogFile, "\n"))

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




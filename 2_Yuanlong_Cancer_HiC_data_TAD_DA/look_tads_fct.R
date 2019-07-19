get_tad_symbols_DT <- function(hicds,
                              exprds,
                              region,
                              mainPipFolder,
                              script0_name = "0_prepGeneData",
                              geneListSuffix = "pipeline_geneList.Rdata", 
                              regionListSuffix = "pipeline_regionList.Rdata",
                              g2tFolder = "genes2tad",
                              g2tSuffix =  "all_genes_positions.txt") {


        
  geneListFile <- file.path(mainPipFolder, hicds, exprds, script0_name, geneListSuffix)
  stopifnot(file.exists(geneListFile))
  geneList <- eval(parse(text = load(geneListFile)))
  
  regionListFile <- file.path(mainPipFolder, hicds, exprds, script0_name, regionListSuffix)
  stopifnot(file.exists(regionListFile))
  regionList <- eval(parse(text = load(regionListFile)))
    
  stopifnot(region %in% regionList)
  
  
  g2tFile <- file.path(hicds, g2tFolder, g2tSuffix)
  stopifnot(file.exists(g2tFile))
  g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
  g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
  
  stopifnot(geneList %in% g2t_DT$entrezID)
  
  stopifnot(region %in% g2t_DT$region)
  
  tad_g2tDT <- g2t_DT[g2t_DT$entrezID %in% geneList &
                      g2t_DT$region == region,
                    ]
  stopifnot(tad_g2tDT$region == region)
  
  tad_entrezID <- as.character(tad_g2tDT$entrezID)
  
  stopifnot(tad_entrezID %in% entrez2symb_dt$entrezID)
  
  tad_symbDT <- entrez2symb_dt[entrez2symb_dt$entrezID %in% tad_entrezID, ]
  
  stopifnot(nrow(tad_symbDT) == length(tad_entrezID))
  
  # sort by name
  tad_symbDT <- tad_symbDT[order(tad_symbDT$symbol),]
  
  outDT <- data.frame(
  hicds = hicds,
  exprds = exprds,
  region = region,
  symbol = tad_symbDT$symbol,
  entrezID = tad_symbDT$entrezID,
  stringsAsFactors = FALSE
  )
  return(outDT)
}

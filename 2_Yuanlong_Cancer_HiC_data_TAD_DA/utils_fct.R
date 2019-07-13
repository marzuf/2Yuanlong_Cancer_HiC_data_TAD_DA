
get_meanCorr_value <- function(exprMatrix, inside_genes, outside_genes, cormet) {
  stopifnot(inside_genes %in% rownames(exprMatrix))
  stopifnot(outside_genes %in% rownames(exprMatrix))
  stopifnot(setequal(c(inside_genes, outside_genes), rownames(exprMatrix)))
  
  nAllGenes <- length(inside_genes) + length(outside_genes)
  
  coexprMatrix <- cor(t(exprMatrix), method = cormet)
  stopifnot(dim(coexprMatrix) == nAllGenes)
  
  coexprMatrix[lower.tri(coexprMatrix, diag = TRUE)] <- NA   # because after I filter that 1 gene should be inside, and 1 gene should be outside -> can never happen to take the diag. value of coexpression
  coexprMatrix <- na.omit(melt(coexprMatrix))
  colnames(coexprMatrix)[1:2] <- c("Var1", "Var2")
  stopifnot(colnames( coexprMatrix)[3] == "value" )
  coexprMatrix$Var1 <- as.character(coexprMatrix$Var1)
  coexprMatrix$Var2 <- as.character(coexprMatrix$Var2)
  
  stopifnot(coexprMatrix$Var1 %in% outside_genes | coexprMatrix$Var1 %in% inside_genes)
  stopifnot(coexprMatrix$Var2 %in% outside_genes | coexprMatrix$Var2 %in% inside_genes)
  stopifnot(inside_genes %in% coexprMatrix$Var1 | inside_genes %in% coexprMatrix$Var2)
  stopifnot(outside_genes %in% coexprMatrix$Var1 | outside_genes %in% coexprMatrix$Var2)
  
  # take only if one of the two genes outside and the other inside
  coexprMatrix <- coexprMatrix[!  (coexprMatrix$Var1 %in% outside_genes & coexprMatrix$Var2 %in% outside_genes),]  # do not take correlation between pairs of genes in  the same TAD
  coexprMatrix <- coexprMatrix[ ! (coexprMatrix$Var1 %in% inside_genes & coexprMatrix$Var2 %in% inside_genes),]  # do not take correlation between pairs of genes in  the same TAD
  
  stopifnot(     (coexprMatrix$Var1 %in% outside_genes & coexprMatrix$Var2 %in% inside_genes) | 
                   (coexprMatrix$Var2 %in% outside_genes & coexprMatrix$Var1 %in% inside_genes) )
  
  

  stopifnot(nrow(coexprMatrix) == (length(inside_genes) * length(outside_genes) ))
  
  meanCorr_value <- mean(coexprMatrix$value)
  stopifnot(!is.na(meanCorr_value))
  return(meanCorr_value)
}

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

checkGenes <- function(genes, side, region_name, gene2tad_dt, tadpos_dt, checkDistKb) {
  
  stopifnot(genes %in% gene2tad_dt$entrezID)
  stopifnot(region_name %in% tadpos_dt$region)
  
  tad_midPos <- (tadpos_dt$start[tadpos_dt$region == region_name]  + tadpos_dt$end[tadpos_dt$region == region_name]  )/2
  stopifnot(length(tad_midPos) == 1)  
    
  genes_midpos <- (gene2tad_dt$start[gene2tad_dt$entrezID %in% genes]+gene2tad_dt$end[gene2tad_dt$entrezID %in% genes])/2
  
  genes_chromo <- gene2tad_dt$chromo[gene2tad_dt$entrezID %in% genes]
  
  stopifnot( genes_chromo == gsub("(chr.+)_TAD.+", "\\1", region_name) )
             
  if(side == "left") {
    stopifnot(genes_midpos <= tad_midPos)
  } else if(side == "right") {
    stopifnot(genes_midpos >= tad_midPos)
  }
    
  if(!is.na(checkDistKb)) {
    
    if(is.numeric(checkDistKb)) {
      checkDistKb_size <- checkDistKb
    } else if(checkDistKb == "sameKb") {
      
      tad_size <- tadpos_dt$end[tadpos_dt$region == region_name]  - tadpos_dt$start[tadpos_dt$region == region_name]  + 1
      stopifnot(length(tad_size) == 1)  
      
      checkDistKb_size <- tad_size
    } else {
      stop("error")
    }
    
    tad_start <- (tadpos_dt$start[tadpos_dt$region == region_name])
    tad_end <- (tadpos_dt$end[tadpos_dt$region == region_name])
    stopifnot(length(tad_start) == 1)  
    stopifnot(length(tad_end) == 1)  
    
    windowLeft <- tad_start - checkDistKb_size
    windowRight <- tad_end + checkDistKb_size 
    
    stopifnot(genes_midpos >= windowLeft)
    stopifnot(genes_midpos <= windowRight)
    
  }
}

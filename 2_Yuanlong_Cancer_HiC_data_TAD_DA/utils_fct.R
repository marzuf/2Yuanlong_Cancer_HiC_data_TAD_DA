source(file.path(setDir, "mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils_meanCorr.R"))

# => get_meanCorr_value() is put in the pipeline folder (TAD_DE_pipeline_v2/TAD_DE_utils_meanCorr.R)
# => get_SAM_FDR_aroundTADs() is put in the pipeline folder (TAD_DE_pipeline_v2/TAD_DE_utils_meanCorr.R)


#############################################################################################################################
############################################################################################################################# 
#############################################################################################################################


get_meanCorr_value_alternative <- function(exprMatrix, sample_genes, cormet) {
  stopifnot(sample_genes %in% rownames(exprMatrix))

  stopifnot(setequal(c(sample_genes), rownames(exprMatrix)))
  
  nAllGenes <- length(sample_genes)
  
  coexprMatrix <- cor(t(exprMatrix), method = cormet)
  stopifnot(dim(coexprMatrix) == nAllGenes)
  
  coexprMatrix[lower.tri(coexprMatrix, diag = TRUE)] <- NA   # because after I filter that 1 gene should be inside, and 1 gene should be outside -> can never happen to take the diag. value of coexpression
  coexprMatrix <- na.omit(melt(coexprMatrix))
  colnames(coexprMatrix)[1:2] <- c("Var1", "Var2")
  stopifnot(colnames( coexprMatrix)[3] == "value" )
  coexprMatrix$Var1 <- as.character(coexprMatrix$Var1)
  coexprMatrix$Var2 <- as.character(coexprMatrix$Var2)
  
  stopifnot(coexprMatrix$Var1 %in% sample_genes & coexprMatrix$Var2 %in% sample_genes)
  

  stopifnot(nrow(coexprMatrix) == (length(sample_genes) * (length(sample_genes)-1) * 0.5 ))
  
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

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

extract_corr_values <- function(sample_values_list, value_to_extract, ds_to_take) {
  stopifnot(ds_to_take %in% names(sample_values_list))
  corr_values <- unlist(lapply(sample_values_list[ds_to_take], function(sub_data) {
    lapply(sub_data, function(x) x[[paste0(value_to_extract)]])
  }))
  corr_values <- na.omit(corr_values)
  return(corr_values)
  
}


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

plot_FDR_with_observedSignif <- function(yaxis_empFDR_vect, xaxis_cutoff, y2_obsSignif, variableName,  feature_name = "genes", rightAxisCol = "steelblue") {
  initmar <- par()$mar
  par(mar = c(5,5,2,5))
  plot((100*yaxis_empFDR_vect) ~ xaxis_cutoff, type="o", pch=16, 
       main=paste0("Empirical FDR with variable cut-off of ", variableName),
       xlab=paste0(variableName, " cut-off"),
       ylab="% FDR", axes=F)
  axis(1, at=xaxis_cutoff, labels = xaxis_cutoff)
  axis(2)
  par(new = T)
  plot(log10(y2_obsSignif) ~ xaxis_cutoff, type="o", 
       col=rightAxisCol, pch=16, ylim=c(0, max(log10(y2_obsSignif))),
       axes=F, xlab=NA, ylab=NA, bty="l")
  box(bty="u")
  axis(side = 4, col=rightAxisCol)
  mtext(side = 4, line = 3, paste0('Observed # of significant ', feature_name, ' (log10)'), col = rightAxisCol)
  on.exit(par(initmar))
}


#############################################################################################################################
############################################################################################################################# 
#############################################################################################################################


# plot_two_FDR_with_observedSignif(first_var_FDR= all_FDR[["empFDR_logFC"]],
#                                    first_var_nbrSignif = all_FDR[["nbrSignif_logFC"]],
#                                              scd_var_FDR = all_FDR[["empFDR_intraTADcorr"]],
#                                  scd_var_nbrSignif = all_FDR[["nbrSignif_intraTADcorr"]],
#                                              first_var_name = "meanFC",
#                                              scd_var_name = "meanCorrV0")
plot_two_FDR_with_observedSignif <- function(first_var_FDR, first_var_nbrSignif,
                                             scd_var_FDR, scd_var_nbrSignif,
                                             first_var_name,
                                             scd_var_name,
                                             first_var_col = "red",
                                             scd_var_col = "blue",
                                             p_pch = 1,
                                             l_lty = 1,
                                             ya_pct = TRUE,
                                             ya_name = "emp. FDR",
                                             ya_name_short = ya_name,
                                             yb_name = "# obs. signif. TADs",
                                             yb_name_short = "# TADs",
                                             yb_log10 = TRUE,
                                             ya_cut_off = NULL,
                                             ya_lty = 2,
                                             ya_col = "darkgrey",
                                             mainLeg = TRUE, 
                                             cutoffLeg = TRUE) {
  stopifnot(names(first_var_FDR) == names(first_var_nbrSignif))
  stopifnot(names(scd_var_FDR) == names(scd_var_nbrSignif)) 
  
  x1 <- as.numeric(names(first_var_FDR))
  x2 <- as.numeric(names(scd_var_FDR))
  stopifnot(!is.na(x1))
  stopifnot(!is.na(x2))
  
  if(yb_log10) {
    y1b <- log10(first_var_nbrSignif)
    y2b <- log10(scd_var_nbrSignif)
    yb_lab <- paste0(yb_name, " (log10)")
  } else {
    y1b <- first_var_nbrSignif
    y2b <- scd_var_nbrSignif
    yb_lab <- paste0(yb_name)
  }
  
  if(ya_pct) {
    y1 <- first_var_FDR*100
    y2 <- scd_var_FDR*100
    y_lab <- paste0(ya_name, " (%)")
    ya_cut_off_line <- ya_cut_off*100
  } else {
    y1 <- first_var_FDR
    y2 <- scd_var_FDR
    y_lab <- paste0(ya_name, " (%)")
    ya_cut_off_line <- ya_cut_off
  }
  y_max <- max(c(y1[!is.infinite(y1)],y2[!is.infinite(y2)]), na.rm = TRUE)
  yb_max <- max(c(y1b[!is.infinite(y1b)],y2b[!is.infinite(y2b)]), na.rm=TRUE)
  
  initmar <- par()$mar
  par(mar = c(5,5,2,5))    
  
  # 1ST VARIABLE: 1st axis
  plot(
    x = x1,
    y = y1,
    type = "l",
    col = first_var_col,
   # pch = 16,
   lty = l_lty,
    axes = FALSE,
    xlab = "",
    ylab = y_lab,
    ylim = c(0, y_max)
  )
  axis(2, at = seq(0, 100, by = 10), labels = TRUE)
  axis(1, line=1,col=first_var_col,col.ticks=first_var_col,col.axis=first_var_col)
  mtext(paste0(first_var_name, " cutoff"),1,line=1,at=mean(x1),col=first_var_col)
  
  if(!is.null(ya_cut_off)) {
    abline(h = ya_cut_off_line, lty = ya_lty, col = ya_col)
  }
  #box(bty="u")
  # 1ST VARIABLE: 2nd axis
  par(new = T)
  plot(
    x = x1,
    y = y1b,
    type = "p",
    col = first_var_col,
    pch = p_pch,
    axes = FALSE,
    xlab = "",
    ylab = paste0(""),
    ylim = c(0, yb_max)
  )
  axis(side = 4, col="black")
  mtext(side = 4, line = 3, paste0(yb_lab), col = "black")
  # 2ND VARIABLE: 1st axis
  par(new = T)
  plot(
    x = x2,
    y = y2,
    type = "l",
    lty = l_lty,
    col = scd_var_col,
   # pch = 16,
    axes = FALSE,
    xlab = "",
    ylab = paste0(""),
    ylim = c(0, y_max)
  )
  axis(1, line=3,col=scd_var_col,col.ticks=scd_var_col,col.axis=scd_var_col)   
  mtext(paste0(scd_var_name, " cutoff"),1,line=3, at=mean(x2),col=scd_var_col)
  # 2ND VARIABLE: 2nd axis
  par(new = T)
  plot(
    x = x2,
    y = y2b,
    type = "p",
    col = scd_var_col,
    pch = p_pch,
    axes = FALSE,
    xlab = "",
    ylab = paste0(""),
    ylim = c(0, yb_max)
  )
  if(mainLeg) {
    legend(x="topright",
           c(ya_name_short, yb_name_short, first_var_name, scd_var_name),
           lty=c(1,-1,-1,-1,-1),
           pch = c(-1, p_pch, -1,-1,-1),
           text.col = c("black", "black", first_var_col, scd_var_col),
           bty="n")
    }
  if(cutoffLeg) {
    legend(x="bottomleft", bty="n", legend=c(paste0(ya_name_short, " cutoff: ", ya_cut_off)), text.col = ya_col, lty = ya_lty, col = ya_col)
  }
  par(mar = initmar)
  on.exit(par(mar = initmar))
}




#############################################################################################################################
############################################################################################################################# 
#############################################################################################################################


plot_meanFC_meanCorr_FDRthresh <- function(
  dataDT,
  var1,
  var2,
  var1_cutoff,
  var2_cutoff,
  var1_name = var1,
  var2_name = var2,
  abs_var1 = FALSE,
  abs_var2 = FALSE,
  annotCol =NULL,
  plotTit = "",
  plotSub = "",
  plotCex = 1.2,
  fileSuffix = NULL,
  plotType = "png",
  legSuppTxt = "",
  mTextLine = -1,
  toPlot = c("grey", "density"),
  fullName=FALSE,
  twoInOne = FALSE
) {

  if(twoInOne) {
   stopifnot(length(toPlot) == 2)
  }

  if(!twoInOne & length(toPlot) == 2)
  stopifnot(!fullName)

  stopifnot(toPlot %in% c("grey", "density"))

  if(is.null(annotCol)) {
    stopifnot(annotCol %in% colnames(dataDT))
  }
  stopifnot(var1 %in% colnames(dataDT))
  stopifnot(var2 %in% colnames(dataDT))
  if(abs_var1) {
    dataDT[, var1] <- abs(dataDT[, var1])
    myxlab <- paste0("abs(", var1_name, ")")
  } else{
    myxlab <- paste0(var1_name)
  }
  if(abs_var2) {
    dataDT[, var2] <- abs(dataDT[, var2])
    myylab <- paste0("abs(", var2_name, ")")
  } else{
    myylab <- paste0(var2_name)
  }
  
  myx <- dataDT[, var1]
  myy <- dataDT[, var2]


 if("density" %in% toPlot) {
  
      ### START PLOTTING THE DENSPLOT
      if(!is.null(fileSuffix)) {
        plotHeight <- ifelse(plotType=="png", 400, 7)
        plotWidth <- plotHeight
        if(twoInOne) plotWidth <- plotWidth*2
        if(fullName) {
          outFile <- paste0(fileSuffix)
        } else {
          if(twoInOne) {
            outFile <- paste0(fileSuffix, "_densAndGreyPlots.", plotType)
          } else {
            outFile <- paste0(fileSuffix, "_densplot.", plotType)
          }  
        }
        dir.create(dirname(outFile), recursive = TRUE)

        do.call(plotType, list(outFile, height = plotHeight, width = plotWidth))
        if(twoInOne) {
          par(mfrow=c(1,2) ) 
        }
      }
      densplot(
        x = myx,
        y = myy,
        xlab = myxlab,
        ylab = myylab,
        cex.lab = plotCex,
        cex.axis = plotCex
      )
      title(main = paste0(plotTit))
      mtext(text = plotSub, side = 3, font=3, line=mTextLine)
      abline(v = var1_cutoff, lty=2, col="darkgrey")
      abline(h = var2_cutoff, lty=2, col="darkgrey")
      legend("bottomright",
             legend = c(legSuppTxt,
                        paste0(var1_name, " cutoff = ", var1_cutoff),
                        paste0(var2_name, " cutoff = ", var2_cutoff)),
             bty="n")
      if(!is.null(fileSuffix)) {
        if(!twoInOne ) {
        foo <- dev.off()
        cat(paste0("... written in densityplot: ", outFile, "\n"))
        }
      }
 }

  if("grey" %in% toPlot) {

      ### START PLOTTING THE DARK GREY PLOT  
      if(!is.null(fileSuffix)) {
        plotHeight <- ifelse(plotType=="png", 400, 7)
        plotWidth <- plotHeight

        if(!twoInOne) {
          if(fullName) {
            outFile <- paste0(fileSuffix)
          } else {
            outFile <- paste0(fileSuffix, "_greyplot.", plotType)
          }
          do.call(plotType, list(outFile, height = plotHeight, width = plotWidth))
        }
      }
      
      if(is.infinite(var1_cutoff) & is.infinite(var2_cutoff)) {
        dotCols <- "black"
      } else if(is.infinite(var2_cutoff)) {
        dotCols <- sapply(1:nrow(dataDT), function(x) ifelse( (dataDT[x,paste0(var1)] >= var1_cutoff), "black", "grey"))
      } else if(is.infinite(var1_cutoff)) {
        dotCols <- sapply(1:nrow(dataDT), function(x) ifelse( (dataDT[x,paste0(var2)] >= var2_cutoff), "black", "grey"))
      } else {
        dotCols <- sapply(1:nrow(dataDT), function(x) ifelse( (dataDT[x,paste0(var1)] >= var1_cutoff & dataDT[x,paste0(var2)] >= var2_cutoff), "black", "grey"))
      }

      plot(
        x = myx,
        y = myy,
        xlab =  myxlab,
        ylab =  myylab,
        pch = 16,
        col = dotCols,
        cex = 0.7,
        cex.lab = plotCex,
        cex.axis = plotCex
      )
      title(main = paste0(plotTit))
      mtext(text = plotSub, side = 3, font=3, line=mTextLine)
      abline(v = var1_cutoff, lty=2, col="red")
      abline(h = var2_cutoff, lty=2, col="red")
      legend("bottomright",
           legend = c(legSuppTxt, 
                      paste0(var1_name, " cutoff = ", var1_cutoff),
                        paste0(var2_name, " cutoff = ", var2_cutoff)),
             bty="n")
      if(!is.null(annotCol)) {
        if(length(dotCols) > 1) {
          toAnnot <- which(dotCols == "black")
          if(length(toAnnot) >= 1)
            text(x = myx[toAnnot],
                 y = myy[toAnnot],
                 labels = dataDT[,paste0(annotCol)][toAnnot],
                 pos=3,
                 cex = 0.7)
        }
      }
      if(!is.null(fileSuffix)) {
        foo <- dev.off()
        cat(paste0("... written in greyplot: ", outFile, "\n"))
      }
  }
}





#############################################################################################################################
############################################################################################################################# 
#############################################################################################################################


plot_meanFC_meanCorr_signifTADs <- function(
  dataDT,
  var1,
  var2,
  signifCol,
  signifTADs,
  var1_name = var1,
  var2_name = var2,
  abs_var1 = FALSE,
  abs_var2 = FALSE,
  annotCol =NULL,
  plotTit = "",
  plotSub = "",
  plotCex = 1.2,
  fileSuffix = NULL,
  plotType = "png",
  legTxt = "",
  mTextLine = -1,
  toPlot = c("grey", "density"),
  fullName=FALSE,
  twoInOne = FALSE
) {

  if(twoInOne) {
   stopifnot(length(toPlot) == 2)
  }

  if(!twoInOne & length(toPlot) == 2)
  stopifnot(!fullName)

  stopifnot(toPlot %in% c("grey", "density"))

  if(is.null(annotCol)) {
    stopifnot(annotCol %in% colnames(dataDT))
  }
  stopifnot(var1 %in% colnames(dataDT))
  stopifnot(var2 %in% colnames(dataDT))
  if(abs_var1) {
    dataDT[, var1] <- abs(dataDT[, var1])
    myxlab <- paste0("abs(", var1_name, ")")
  } else{
    myxlab <- paste0(var1_name)
  }
  if(abs_var2) {
    dataDT[, var2] <- abs(dataDT[, var2])
    myylab <- paste0("abs(", var2_name, ")")
  } else{
    myylab <- paste0(var2_name)
  }
  
  myx <- dataDT[, var1]
  myy <- dataDT[, var2]


 if("density" %in% toPlot) {
  
      ### START PLOTTING THE DENSPLOT
      if(!is.null(fileSuffix)) {
        plotHeight <- ifelse(plotType=="png", 400, 7)
        plotWidth <- plotHeight
        if(twoInOne) plotWidth <- plotWidth*2
        if(fullName) {
          outFile <- paste0(fileSuffix)
        } else {
          if(twoInOne) {
            outFile <- paste0(fileSuffix, "_densAndGreyPlots.", plotType)
          } else {
            outFile <- paste0(fileSuffix, "_densplot.", plotType)
          }  
        }
        dir.create(dirname(outFile), recursive = TRUE)
        cat("... open in density\n")
        do.call(plotType, list(outFile, height = plotHeight, width = plotWidth))
        if(twoInOne) {
          par(mfrow=c(1,2) ) 
        }
      }
      densplot(
        x = myx,
        y = myy,
        xlab = myxlab,
        ylab = myylab,
        cex.lab = plotCex,
        cex.axis = plotCex
      )
      title(main = paste0(plotTit))
      mtext(text = plotSub, side = 3, font=3, line=mTextLine)

      legend("bottomright",
           legend = c(legTxt),
             bty="n")


      if(!is.null(fileSuffix)) {
        if(!twoInOne ) {
        foo <- dev.off()
        cat(paste0("... written densityplot: ", outFile, "\n"))
        }
      }
 }

  if("grey" %in% toPlot) {

      ### START PLOTTING THE DARK GREY PLOT  
      if(!is.null(fileSuffix)) {
        plotHeight <- ifelse(plotType=="png", 400, 7)
        plotWidth <- plotHeight

        if(!twoInOne) {
          if(fullName) {
            outFile <- paste0(fileSuffix)
          } else {
            outFile <- paste0(fileSuffix, "_greyplot.", plotType)
          }
          cat("... open in grey\n")
          do.call(plotType, list(outFile, height = plotHeight, width = plotWidth))
        }
      }
      


      dotCols <- ifelse(dataDT[,signifCol] %in% signifTADs, "black", "grey")

      plot(
        x = myx,
        y = myy,
        xlab =  myxlab,
        ylab =  myylab,
        pch = 16,
        col = dotCols,
        cex = 0.7,
        cex.lab = plotCex,
        cex.axis = plotCex
      )
      title(main = paste0(plotTit))
      mtext(text = plotSub, side = 3, font=3, line=mTextLine)

      legend("bottomright",
           legend = c(legTxt),
             bty="n")

      if(!is.null(annotCol)) {
        if(length(dotCols) > 1) {
          toAnnot <- which(dotCols == "black")
          if(length(toAnnot) >= 1)
            text(x = myx[toAnnot],
                 y = myy[toAnnot],
                 labels = dataDT[,paste0(annotCol)][toAnnot],
                 pos=3,
                 cex = 0.7)
        }
      }
      if(!is.null(fileSuffix)) {
        foo <- dev.off()
        cat(paste0("... written greyplot: ", outFile, "\n"))
      }
  }
}








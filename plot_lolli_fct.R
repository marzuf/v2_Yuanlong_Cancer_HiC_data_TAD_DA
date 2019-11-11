
######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################

plot_lolliTAD <- function(TAD_to_plot, 
                          meanExpr_vect,
                          DE_table,
                          g2t_table,
                          id2name_table,
                          geneList,
						  barcolBy = "logFC", 
                          orderBy = "startPos",
                          upcolor = "limegreen", 
                          downcolor = "orangered",
                          palettefunc_pos = colorRampPalette(c("#e6ffe6", "#00b300")),
                          palettefunc_neg = colorRampPalette(c("#ffe6e6", "#b30000")), 
                          graybars = FALSE,
                          textLeft = FALSE,
						  cond1 = "cond1",
						  cond2 = "cond2",
						  labelWithRank=FALSE,
						  mytitle=NULL
                          ) {

  suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  suppressPackageStartupMessages(library(plotrix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  suppressPackageStartupMessages(library(gridExtra, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

  stopifnot(barcolBy %in% c("meanExpr", "logFC"))


library(rlang)
as_dictionary <- as_data_pronoun

  stopifnot("entrezID" %in% colnames(id2name_table))
  stopifnot("geneName" %in% colnames(id2name_table))
  
  rownames(id2name_table) <- id2name_table$entrezID
  rownames(g2t_table) <- g2t_table$entrezID
  rownames(DE_table) <- DE_table$genes



  if(!orderBy %in% c("logFC", "meanExpr_vect", "startPos"))
    stop("ERROR - invalid \"orderBy\" argument !\n")
  negPalette <- palettefunc_neg(11)
  posPalette <- palettefunc_pos(11)
  
  genes_to_plot <- g2t_table$entrezID[g2t_table$region == TAD_to_plot]
  genes_to_plot <- genes_to_plot[genes_to_plot %in% geneList]
  if(length(genes_to_plot) == 0) {
    warning("ERROR - no genes to plot \n")
	return(grob())
   }
  TAD_genes_DT <- data.frame(gene = genes_to_plot, 
                             symbol = id2name_table[genes_to_plot, "geneName"],
                             start = g2t_table[genes_to_plot, "start"],
                             mean_expr = meanExpr_vect[genes_to_plot],
                             log_FC = DE_table[genes_to_plot, "logFC"],
                             voom_adj_pval = DE_table[genes_to_plot, "adj.P.Val"])
  
  stopifnot(!any(is.na(TAD_genes_DT)))
  
  TAD_genes_DT$start <- as.numeric(as.character(TAD_genes_DT$start))
  TAD_genes_DT$gene <- as.character(TAD_genes_DT$gene)
  
  if(orderBy == "logFC") {
    TAD_genes_DT <- TAD_genes_DT[order(TAD_genes_DT$log_FC),]
  } else if(orderBy == "meanExpr") {
    TAD_genes_DT <- TAD_genes_DT[order(TAD_genes_DT$mean_expr),]
  } else if(orderBy == "startPos") {
    TAD_genes_DT <- TAD_genes_DT[order(TAD_genes_DT$start),]
  } else {
    stop("error")
  }   
  
  TAD_genes_DT$gene <- factor(TAD_genes_DT$gene, levels = as.character(TAD_genes_DT$gene))
  
  TAD_genes_DT$signif <-  unlist(sapply(TAD_genes_DT$voom_adj_pval, function(x)
        ifelse(x < 0.001, "***", 
               ifelse(x < 0.01, "**", 
                      ifelse(x < 0.05, "*", "")))))
  
  TAD_genes_DT$logFC_color <- unlist(sapply(TAD_genes_DT$log_FC, function(x){
    if(x < 0)
      return(downcolor)
    if(x > 0) 
      return(upcolor)
    return("black")
  }))
  
  ### ADDED FOR THE RANK - 11.04.18
  DE_table <- DE_table[order(DE_table$adj.P.Val),]

  TAD_genes_DT$gene_rank <-  unlist(sapply(TAD_genes_DT$gene, function(x) which(as.character(DE_table$genes) == as.character(x)) ))



  
  if(graybars) {
    barColor <- "lightgray"
  } else{
	if(barcolBy == "meanExpr") {
	  resc_mean_expr <- plotrix::rescale(TAD_genes_DT$mean_expr, newrange=c(0,1))
      TAD_genes_DT$meanExpr_color <- unlist(foreach(i = 1:nrow(TAD_genes_DT), .combine='c') %do% {
        exprVal <- resc_mean_expr[i]
	    ifelse(TAD_genes_DT$log_FC[i] < 0, negPalette[round(exprVal*10)+1], 
    	       ifelse(TAD_genes_DT$log_FC[i] > 0, posPalette[round(exprVal*10)+1], "gray"))
		  })
	  barColor <- TAD_genes_DT$meanExpr_color
	} else if(barcolBy == "logFC") {

	  resc_logFC <- plotrix::rescale(abs(TAD_genes_DT$log_FC), newrange=c(0,1))
      TAD_genes_DT$logFC_barcolor <- unlist(foreach(i = 1:nrow(TAD_genes_DT), .combine='c') %do% {
        logVal <- resc_logFC[i]
	    ifelse(TAD_genes_DT$log_FC[i] < 0, negPalette[round(logVal*10)+1], 
    	       ifelse(TAD_genes_DT$log_FC[i] > 0, posPalette[round(logVal*10)+1], "gray"))
		  })
	  barColor <- TAD_genes_DT$logFC_barcolor
	} else {
	  stop("---error\n")
	}
  }
  
  my_xlab <- ifelse(labelWithRank, paste0("# genes =  ",  nrow(DE_table)), "")

#cat("AAA\n")
  
  lolliTitle <- ifelse(is.null(mytitle), TAD_to_plot, mytitle) 
  
  
  # use gene entrez ID not symbols here (some genes duplicated name with different entrez ID)
  p <- ggdotchart(TAD_genes_DT, x = "gene", y = "log_FC",
                  title = lolliTitle,
             color = TAD_genes_DT$logFC_color ,             # color of the dots
             ylab = "Log2(fold change)",
             xlab=my_xlab,
              add = "segments",                             # Add segments from y = 0 to dots
             add.params = list(color = barColor, size = 2), # Change segment color and size
             
             # for correct sorting
             order = as.character(TAD_genes_DT$gene),
             sort.by.groups = FALSE,
             group = "gene",
             dot.size = 10,                                 # Large dot size
             label = TAD_genes_DT$signif,                        # Add mpg values as dot labels
              font.label = list(color = "white", size = 14, 
                               vjust = 0.5),               # Adjust label parameters
             ggtheme = theme_pubr()                        # ggplot2 theme
  ) + 
    geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
    theme(plot.title = element_text(hjust = 0.5, face=2, size=18),
          axis.title.x = element_text(hjust = 0.5, vjust = 1, face=3, size=10)
    )

#cat("BBB\n")
  
# change the x labels to be the gene names (ok because TAD_genes_DT has been sorted !!!)
 if(labelWithRank){
   TAD_genes_DT$symbol <- paste0(TAD_genes_DT$symbol, "\n(", TAD_genes_DT$gene_rank, ")")
 }
 p <- p+scale_x_discrete(labels=TAD_genes_DT$symbol)
  
  p <- ggpar(p, ylim=c(-max(abs(TAD_genes_DT$log_FC), na.rm=T), max(abs(TAD_genes_DT$log_FC), na.rm=T)))
  p <- p + annotate("text", x = ifelse(textLeft, 0 , nrow(TAD_genes_DT))+0.5, y = max(abs(TAD_genes_DT$log_FC), na.rm=T), 
                    label = paste0("expr. ", cond1, " >\nexpr. ", cond2), color=upcolor, hjust =  ifelse(textLeft,0,1), fontface="bold")
  p <- p + annotate("text", x = ifelse(textLeft, 0 , nrow(TAD_genes_DT))+0.5, y = -max(abs(TAD_genes_DT$log_FC), na.rm=T), 
                    label = paste0("expr. ", cond2, " >\nexpr. ", cond1), color=downcolor, hjust = ifelse(textLeft,0,1), fontface="bold")
  return(p)
}



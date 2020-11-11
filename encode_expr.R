# Rscript encode_expr.R

require(foreach)
require(doMC)
registerDoMC(40)

require(ggpubr)

dataFolder <- "encode_rnaseq"
# dataFile <- "ENCSR444WCZ_A549_ENCFF133PXI.tsv"
# dataFile <- "ENCSR489OCU_NCI-H460_ENCFF562TXM.tsv"

outFolder <- "ENCODE_EXPR"
dir.create(outFolder, recursive = TRUE)

id_of_interest <- c("1645", "1646") # AKR1C1 AKR1C2
ids_of_interest <- c("AKR1C1" = "ENSG00000187134", "AKR1C2"= "ENSG00000151632", "AKR1C3"="ENSG00000196139") # AKR1C1 AKR1C2
ensg2symb<- c( "ENSG00000187134" ="AKR1C1", "ENSG00000151632" = "AKR1C2", "ENSG00000196139"="AKR1C3") # AKR1C1 AKR1C2


cl_names <- c("ENCSR489OCU_NCI-H460" = "H460", "ENCSR444WCZ_A549" = "A549", "LG1" = "LG1", "LG2"="LG2")

count_col <- "FPKM"
# count_col <- "TPM"
# count_col <- "expected_count"

aggFun <- "mean"

buildTable <- TRUE


all_files <- list.files(dataFolder, pattern="tsv$")

if(buildTable) {
  
  all_plot_dt <- foreach(dataFile = all_files, .combine='rbind') %dopar% {
    
    cl <- gsub("(.+)_(EN.+tsv)", "\\1", dataFile)
    
    
    in_dt <- read.delim(file.path(dataFolder, dataFile), stringsAsFactors = FALSE, header=T)
    
    
    
    to_keep <- sapply(in_dt$gene_id, function(x) any(sapply(ids_of_interest, function(genepatt) grepl(genepatt, x))))
    
    oi_dt <- in_dt[which(to_keep),]
    
    cat(paste0("... nrow before agg: \t", nrow(oi_dt), "\n"))
    
    stopifnot(nrow(oi_dt) == sum(sapply(ids_of_interest, function(x) sum(grepl(x, in_dt$gene_id)))))
    
    
    oi_dt$gene_id_short <- gsub("(ENSG.+)\\..+$", "\\1", oi_dt$gene_id)
    
    stopifnot(oi_dt$gene_id_short %in%ids_of_interest)
    
    
    plot_dt <- aggregate(as.formula(paste0(count_col, "~gene_id_short")), data = oi_dt, FUN=aggFun)
    
    plot_dt$cl <- cl
    plot_dt
    
  }
  
  outFile <- file.path(outFolder, "all_plot_dt.Rdata")
  save(all_plot_dt, file=outFile, version=2)
  cat(paste0("... written:", outFile,"\n"))
  
  
  
} else {
  outFile <- file.path(outFolder, "all_plot_dt.Rdata")
  all_plot_dt <- get(load(outFile))
  # load("ENCODE_EXPR/all_plot_dt.Rdata")
  
  
  
}


all_plot_dt$symbol <- ensg2symb[all_plot_dt$gene_id_short]
stopifnot(!is.na(all_plot_dt$symbol))
all_plot_dt$cl_lab <- cl_names[all_plot_dt$cl]
stopifnot(!is.na(all_plot_dt$cl_lab))


mygene =unique(all_plot_dt$symbol)[1]

for(mygene in unique(all_plot_dt$symbol)){
  
  sub_dt <- all_plot_dt[all_plot_dt$symbol == mygene,]

  p <- ggboxplot(sub_dt, y = paste0(count_col), x ="cl_lab", add="jitter",
            xlab ="") + 
    ggtitle(paste0(mygene), subtitle=paste0(count_col)) +
    # scale_color_manual(values="my_cols")+
    # scale_fill_manual(values=my_cols)  +
    # labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
    # guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
      # text = element_text(family=fontFamily),
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
  
  outFile <- file.path(outFolder, paste0(mygene, "_", count_col, "_encode_geneeexpr_boxplot.", "svg"))
  ggsave(p, file=outFile, height=5.5, width=7)
  cat(paste0("... written: ", outFile, "\n"))
  
    
}





setDir <- "/media/electron"
setDir <- ""

# Rscript revision_top4_expr_cell_lines.R
script_name="revision_top4_expr_cell_lines.R"
startTime <- Sys.time()
cat("> START ", script_name, "\n")

hicds1 <- "GSE118514_RWPE1_40kb"
exprds <- "TCGAprad_norm_prad"
all_tads <- "chr12_TAD194"
all_tads <- "chr17_TAD147"
hicds2 <- "GSE118514_22Rv1_40kb"

plotType <- "svg"
myHeightGG <- 5
myWidthGG <- 6

# Rscript revision_top4_expr_cell_lines.R GSE118514_RWPE1_40kb TCGAprad_norm_prad GSE118514_22Rv1_40kb chr12_TAD194 chr17_TAD147 chr7_TAD424 chr17_TAD174

# Rscript revision_top4_expr_cell_lines.R GSE118514_22Rv1_40kb TCGAprad_norm_prad GSE118514_RWPE1_40kb chr12_TAD196 chr1_TAD460 chr17_TAD135 chr17_TAD269

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 3) {
  hicds1 <- args[1]
  exprds <- args[2]
  hicds2 <- args[3]
  all_tads <- args[4:length(args)]
}

cat(paste0("hicds1 = ", hicds1, "\n"))
cat(paste0("hicds2 = ", hicds2, "\n"))

require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
# minGenes <- 3
expr_var <- "FPKM"

log10_offset <- 0.01
plot_var <- paste0(expr_var, "_log10")

outFolder <- "REVISION_TOP4_EXPR_CELL_LINES"
dir.create(outFolder, recursive = TRUE)

#### will not look at mean corr because max number of replicates is 6, with median=2

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

# symb_aggFun <- "mean"  # aggregate expression data to have unique gene symbol
# rep_aggFun <- "mean"  # aggregate expression across replicates
# tad_aggFun <- "mean" # aggregate expression to TAD level

# all_mrna_data=get(load("sub_expression_allDS.Rdata"))

all_mrna_data <-get(load(file.path(setDir,
                                   "/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/1.Data/exp_CELL_LINEs/exp_info.Rdata")))

# all_mrna_data <- get(load("sub_all_mrna_data.Rdata"))

source("revision_settings.R")

# all_hicds=all_hicds[1]


stopifnot(hicds1 %in% names(all_hicds))
ds1 <- as.character(all_hicds[paste0(hicds1)])
stopifnot(length(ds1) > 0)
stopifnot(!is.na(ds1))
stopifnot(ds1 %in% names(all_mrna_data))
ds1_mrna_data <- all_mrna_data[[ds1]]


stopifnot(hicds2 %in% names(all_hicds))
ds2 <- as.character(all_hicds[paste0(hicds2)])
stopifnot(length(ds2) > 0)
stopifnot(!is.na(ds2))
stopifnot(ds2 %in% names(all_mrna_data))
ds2_mrna_data <- all_mrna_data[[ds2]]



## laod the gene2tad assignment
region_file <- file.path(hicds1, "genes2tad", "all_genes_positions.txt")
genes2tad_dt <- read.delim(region_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID",  "chromo", "start", "end", "region"))
genes2tad_dt$entrezID <- as.character(genes2tad_dt$entrezID)                                                               
stopifnot(genes2tad_dt$entrezID %in% names(entrez2symb))
genes2tad_dt$symbol <- entrez2symb[paste0(genes2tad_dt$entrezID)]
stopifnot(!is.na(genes2tad_dt$symbol))
# take only TAD regions
genes2tad_dt <- genes2tad_dt[grepl("_TAD", genes2tad_dt$region),]

# take only the genes included in the pipeline
pipeline_geneList <- get(load(file.path("PIPELINE/OUTPUT_FOLDER/", hicds1, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
stopifnot(pipeline_geneList %in% genes2tad_dt$entrezID)
genes2tad_dt <- genes2tad_dt[genes2tad_dt$entrezID %in% pipeline_geneList,]

for(tad_to_plot in all_tads) {
  
  cat(paste0("for ", hicds1, " - ", tad_to_plot, " (vs ", hicds2, ")"))
  
  stopifnot(tad_to_plot %in% genes2tad_dt$region)
  tad_symbols <- genes2tad_dt$symbol[genes2tad_dt$region == tad_to_plot]
  stopifnot(length(tad_symbols) > 0)
  
  if("PRAC" %in% tad_symbols) tad_symbols <- c("PRAC1", "PRAC2", tad_symbols)
  
  # tad_symbols <- c("FGR", "CFH")
  
  # iterate over replicates for the 1st hicds
  
  hicds1_expr_dt <- do.call(rbind, lapply(1:length(ds1_mrna_data), function(i_x) {
    x <- ds1_mrna_data[[i_x]]
    
    stopifnot(any(x$external_gene_name %in% genes2tad_dt$symbol))
    
    stopifnot("external_gene_name" %in% colnames(x))
    stopifnot(paste0(expr_var) %in% colnames(x))
    sub_x <- x[x$external_gene_name %in% tad_symbols,,drop=F]
    if(nrow(sub_x) == 0) {
      sub_x$rep <- character(0)
      sub_x$hicds <- character(0)
    } else {
      sub_x$rep <- names(ds1_mrna_data)[i_x]
      sub_x$hicds <- hicds1
    }
    sub_x
  }))
  cat(paste0("nrow hicds1 = ", nrow(hicds1_expr_dt), "\n"))
  
  
  hicds2_expr_dt <- do.call(rbind, lapply(1:length(ds2_mrna_data), function(i_x) {
    x <- ds2_mrna_data[[i_x]]
    stopifnot("external_gene_name" %in% colnames(x))
    stopifnot(paste0(expr_var) %in% colnames(x))
    sub_x <- x[x$external_gene_name %in% tad_symbols,,drop=F]
    if(nrow(sub_x) == 0) {
      sub_x$rep <- character(0)
      sub_x$hicds <- character(0)
    } else {
      sub_x$rep <- names(ds1_mrna_data)[i_x]
      sub_x$hicds <- hicds2
    }
    sub_x
  }))
  
  cat(paste0("nrow hicds2 = ", nrow(hicds2_expr_dt), "\n"))
  
  outFile <- file.path(outFolder, paste0(hicds1, "_", tad_to_plot, "_hicds1_expr_dt.Rdata"))
  save(hicds1_expr_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, paste0(hicds1, "_", tad_to_plot, "_hicds2_expr_dt.Rdata"))
  save(hicds2_expr_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  keepcols <- c("hicds", "rep", expr_var, "external_gene_name")
  tad_expr_dt <- rbind(hicds1_expr_dt[,keepcols], hicds2_expr_dt[,keepcols])
  
  tad_expr_dt[,paste0(expr_var, "_log10")] <- log10(tad_expr_dt[,paste0(expr_var)]+log10_offset)
  
  
  tad_expr_dt$hicds <- gsub("_40kb", "", tad_expr_dt$hicds)
  plot_dt <- tad_expr_dt
  outFile <- file.path(outFolder, paste0(hicds1, "_", tad_to_plot, "_", hicds2, "_", plot_var, "_dt.Rdata"))
  save(plot_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  mycols <- c("steelblue3", "tan2")
  my_xlab <- ""
  my_ylab <- paste0("mRNA ", expr_var, "(log10[+", log10_offset, "])")
  plotTit <- paste0(hicds1, " - " , tad_to_plot, ": mRNA (vEntrez)")
  subTit <- paste0("# rep. ", gsub("_40kb", "",hicds1), " = ", length(unique(hicds1_expr_dt$rep)), "; # rep. ",gsub("_40kb", "",hicds2), " = ",  length(unique(hicds2_expr_dt$rep)))
  
  p_var_boxplot <- ggplot(tad_expr_dt, aes_string(x = "external_gene_name", y = plot_var, color = "hicds")) + 
    geom_boxplot(notch = F, outlier.shape=NA)+
    geom_point(aes(fill=hicds), position=position_jitterdodge(),  alpha=0.5) +
    
    ggtitle(paste0(plotTit), subtitle = paste0(subTit))+
    scale_x_discrete(name=my_xlab)+
    scale_y_continuous(name=paste0(my_ylab),
                       breaks = scales::pretty_breaks(n = 20))+
    
    scale_color_manual(values=mycols)+
    scale_fill_manual(values=mycols)+
    
    labs(fill  = paste0("DS"), color=paste0("DS")) +
    theme( 
      # text = element_text(family=fontFamily),
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      panel.grid.minor.y = element_line(colour = "grey"),
      axis.line.x= element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .2, color = "black"),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
      axis.text.x =element_text(color="black", hjust=0.5,vjust = 0.5, size=12, face="bold"),
      # axis.ticks.x = element_blank(),
      axis.title.y = element_text(color="black", size=13),
      axis.title.x = element_text(color="black", size=13),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.text = element_text(size=12),
      legend.key = element_blank(),
      legend.title = element_text(face="bold")
    )
  
  outFile <- file.path(outFolder, paste0(hicds1, "_", tad_to_plot, "_", hicds2, "_", plot_var, "_boxplot.", plotType))
  ggsave(plot = p_var_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
  cat(paste0("... written: ", outFile, "\n"))
  
  
}

  
# 
# tad_expr_dt <- merge(hicds1_expr_dt[,keepcols],
#                       hicds2_expr_dt[,keepcols],
#                       by=c("external_gene_name"), all=FALSE)

# subTit <- ""
# 
# 
# ggboxplot(tad_expr_dt, 
#                x = "hicds",
#                xlab = paste0(""), 
#                y = paste0(plot_var),
#                ylab = paste0("mRNA ", expr_var, "(log10[+", log10_offset, "])"),
#                legend.title="",
#                # legend="right",
#                legend="none",
#                title = paste0(tad_to_plot, ": mRNA (vEntrez)"),
#                subtitle=subTit,
#                fill = "hicds", palette = c("steelblue3", "tan2"),
#                add = "jitter") + 
#   facet_wrap(~external_gene_name)
# 
# 
# p <- p + font("legend.text", size = 14)
# p <- p + font("subtitle", face="bold")
# p <- p + font("title", size = 18, face="bold")
# p <- p + font("xlab", size = 16, face="bold")
# p <- p + font("ylab", size = 16)
# p <- p + font("xy.text", size = 14)
# p <- p + theme(plot.title = element_text(hjust=0.5))
# p <- p + stat_compare_means(label.x=0.5)
# if(SSHFS) p
# 
# outFile <- file.path(outFold, paste0(curr_dataset, "_", exprds, "_", tad, "_mRNA", countType, "_vEntrez.", plotType))
# ggsave(plot=p, filename = outFile, height=myHeight, width=myWidth )
# cat(paste0("... written: ", outFile, "\n"))

######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

# tmp <- list(c("A"=1,"B"=2,"C"=3), c("A"=3,"B"=6,"C"=9), c("A"=2, "B"=1, "C" =3))
# Reduce(rep_aggFun, tmp)  ### !!! not working
# Reduce(`+`, tmp)/length(tmp)
# colMeans(do.call(rbind, tmp))
# all(colMeans(do.call(rbind, tmp)) == Reduce(`+`, tmp)/length(tmp))
#
# apply( do.call(rbind, tmp), 2,rep_aggFun)
# all(apply( do.call(rbind, tmp), 2,rep_aggFun) == Reduce(`+`, tmp)/length(tmp))

# work only if length >= 2
# all(sapply(lapply(expr_values[2:length(expr_values)], names), FUN = identical, names(expr_values[[1]])))
# all(sapply(lapply(tmp[2:length(tmp)], names), FUN = identical, names(tmp[[1]])))
# > tmp <- list(c("A"=1,"B"=2,"C"=3), c("A"=3,"B"=6,"C"=9))
# > unique(lapply(tmp, names))
# [[1]]
# [1] "A" "B" "C"
#
# > Reduce(`+`, tmp)
# A  B  C
# 4  8 12

# > tmp <- list(c(1,2,3), c(1,2,3))
# > unique(tmp)
# [[1]]
# [1] 1 2 3
#
# > tmp <- list(c(1,2,3), c(1,2,4))
# > unique(tmp)
# [[1]]
# [1] 1 2 3
#
# [[2]]
# [1] 1 2 4
#
# >

# > sub_data2=x[["mega_GSE105318_DLD1"]]
# > str(sub_data2)
# List of 2
# $ rep1:Classes ‘data.table’ and 'data.frame':  54849 obs. of  9 variables:
#   ..$ ensembl_gene_id   : chr [1:54849] "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" "ENSG00000000938" ...
# ..$ TPM               : num [1:54849] 79.87 2.28 17.89 0 0.17 ...
# ..$ FPKM              : num [1:54849] 70.91 2.03 15.88 0 0.15 ...
# ..$ expected_count    : num [1:54849] 499 81.7 256.3 0 2 ...
# ..$ external_gene_name: chr [1:54849] "DPM1" "SCYL3" "C1orf112" "FGR" ...
# ..$ chromosome_name   : num [1:54849] 20 1 1 1 1 6 6 6 1 1 ...
# ..$ start_position    : int [1:54849] 49551404 169818772 169631245 27938575 196621008 143815948 53362139 41040684 24683489 24742284 ...
# ..$ end_position      : int [1:54849] 49575092 169863408 169823221 27961788 196716634 143832827 53481768 41067715 24743424 24799466 ...
# ..$ gene_biotype      : chr [1:54849] "protein_coding" "protein_coding" "protein_coding" "protein_coding" ...
# ..- attr(*, ".internal.selfref")=<externalptr>
#   ..- attr(*, "sorted")= chr "ensembl_gene_id"
# $ rep2:Classes ‘data.table’ and 'data.frame':  54849 obs. of  9 variables:
#   ..$ ensembl_gene_id   : chr [1:54849] "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" "ENSG00000000938" ...
# ..$ TPM               : num [1:54849] 88.64 2.61 16.89 0 0.2 ...
# ..$ FPKM              : num [1:54849] 75.12 2.21 14.32 0 0.17 ...
# ..$ expected_count    : num [1:54849] 627 105 279 0 2 ...
# ..$ external_gene_name: chr [1:54849] "DPM1" "SCYL3" "C1orf112" "FGR" ...
# ..$ chromosome_name   : num [1:54849] 20 1 1 1 1 6 6 6 1 1 ...
# ..$ start_position    : int [1:54849] 49551404 169818772 169631245 27938575 196621008 143815948 53362139 41040684 24683489 24742284 ...
# ..$ end_position      : int [1:54849] 49575092 169863408 169823221 27961788 196716634 143832827 53481768 41067715 24743424 24799466 ...
# ..$ gene_biotype      : chr [1:54849] "protein_coding" "protein_coding" "protein_coding" "protein_coding" ...
# ..- attr(*, ".internal.selfref")=<externalptr>
#   

# -> for each TAD
# -> how many genes available in the  rna expr
# -> mean FPKM
# -> then I can do the same as for norm and tumor
# Danielle prefer TPM, but I use FPKM


# Rscript pands_coreg_selectedTADs.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD390

outFolder <- "PANDS_COREG_SELECTEDTADS"
dir.create(outFolder, recursive = TRUE)

hicds_oi <- "ENCSR489OCU_NCI-H460_40kb"
exprds_oi <- "TCGAluad_norm_luad"
tad_oi <- "chr11_TAD390"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
hicds_oi <- args[1]
exprds_oi <- args[2]
tad_oi <- args[3]

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"

setDir <- "/media/electron"
setDir <- ""
buildTable <- FALSE

corrMeth <- "pearson"

plotType <- "svg"
myWidth <- myHeight <- 8

require(doMC)
require(foreach)
registerDoMC(40)
require(ggplot2)

aggFun <- "mean" # aggregate correlation for the other datasets


pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds


### retrieve the genes that belong to the TAD of interest
pipeline_geneList <- get(load(file.path(pipFolder, hicds_oi, exprds_oi, script0_name, "pipeline_geneList.Rdata")))
gene2tad_dt <- read.delim(file.path(hicds_oi, "genes2tad", "all_genes_positions.txt"),
                          stringsAsFactors = FALSE,
                          header=FALSE, 
                          col.names=c("entrezID", "chromosome", "start", "end", "region"))
gene2tad_dt$entrezID <- as.character(gene2tad_dt$entrezID)
stopifnot(pipeline_geneList %in% gene2tad_dt$entrezID)
stopifnot(tad_oi %in% gene2tad_dt$region)
g2t_dt <- gene2tad_dt[gene2tad_dt$entrezID %in% pipeline_geneList,]
genes_oi <- g2t_dt$entrezID[g2t_dt$region == tad_oi]
all_gene_pairs <- combn(genes_oi, m=2)

rm("g2t_dt")
rm("gene2tad_dt")
rm("pipeline_geneList")

### retrieve pairwise correlations for all datasets (only if they are in the same TAD)

if(buildTable) {
  hicds = all_hicds[1]
  all_ds_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat(paste0("... start: ", hicds, " - ", exprds, "\n"))
      source(file.path("PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R")))
      samp1 <- get(load(file.path(setDir, sample1_file)))
      samp2 <- get(load(file.path(setDir, sample2_file)))
      qq_dt <- eval(parse(text = load(file.path(pipFolder, hicds, exprds, script0_name, "rna_qqnorm_rnaseqDT.Rdata")))) 
      stopifnot(setequal(colnames(qq_dt), c(samp1,samp2)))
      pipeline_geneList <- get(load(file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")))
      gene2tad_dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"),
                                stringsAsFactors = FALSE,
                                header=FALSE, 
                                col.names=c("entrezID", "chromosome", "start", "end", "region"))
      gene2tad_dt$entrezID <- as.character(gene2tad_dt$entrezID)
      stopifnot(pipeline_geneList %in% gene2tad_dt$entrezID)
      g2t_dt <- gene2tad_dt[gene2tad_dt$entrezID %in% pipeline_geneList,]
      do.call(rbind, apply(all_gene_pairs, 2, function(gp) {
        gene1 <- gp[1]
        gene2 <- gp[2]
        tad_gene1 <- g2t_dt$region[g2t_dt$entrezID == gene1]
        tad_gene2 <- g2t_dt$region[g2t_dt$entrezID == gene2]
        if(length(tad_gene1) == 0 | length(tad_gene2) == 0) {
          pcorr <- NA
        }else {
          if(!gene1 %in% rownames(qq_dt) | !gene2 %in% rownames(qq_dt)) {
            pcorr <- NA
          } else {
            if(tad_gene1 == tad_gene2) {
              pcorr <- cor(qq_dt[gene1,], qq_dt[gene2,], method=corrMeth)  
            } else {
              pcorr <- NA
            }
          }
        }
        data.frame(
          hicds=hicds,
          exprds=exprds,
          gene1=gene1,
          gene2=gene2,
          corr=pcorr,
          stringsAsFactors = FALSE
        )
      }))
      
    }
    exprds_dt
  }
  
  outfile <-file.path(outFolder, paste0(hicds_oi, "_", exprds_oi, "_", tad_oi, "_all_ds_dt.Rdata"))
  save(all_ds_dt, file=outfile, version=2)
  cat(paste0("... written:" , outfile, "\n"))
  
} else {
  outfile <-file.path(outFolder, paste0(hicds_oi, "_", exprds_oi, "_", tad_oi, "_all_ds_dt.Rdata"))
  all_ds_dt <- get(load(outfile))
}
# load("PANDS_COREG_SELECTEDTADS/all_ds_dt.Rdata")
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(is.numeric(gff_dt$end))
stopifnot(is.numeric(gff_dt$start))
gff_dt <- gff_dt[order(gff_dt$chromo, gff_dt$start, gff_dt$end),]
symbols_level <- gff_dt$symbol[gff_dt$entrezID %in% genes_oi]

stopifnot(all_ds_dt$gene1 %in% gff_dt$entrezID)
stopifnot(all_ds_dt$gene2 %in% gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$symbol))
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol,gff_dt$entrezID)
all_ds_dt$symbol1 <- entrez2symb[paste0(all_ds_dt$gene1)]
all_ds_dt$symbol2 <- entrez2symb[paste0(all_ds_dt$gene2)]
stopifnot(!is.na(all_ds_dt$symbol1))
stopifnot(!is.na(all_ds_dt$symbol2))

oi_ds_dt <- all_ds_dt[all_ds_dt$hicds == hicds_oi & all_ds_dt$exprds == exprds_oi,]
stopifnot(nrow(oi_ds_dt) == length(genes_oi) * (length(genes_oi)-1) * 0.5 )
stopifnot(!is.na(oi_ds_dt))

other_ds_dt <- all_ds_dt[!(all_ds_dt$hicds == hicds_oi & all_ds_dt$exprds == exprds_oi),]
stopifnot(nrow(oi_ds_dt) + nrow(other_ds_dt) == nrow(all_ds_dt))
other_ds_dt <- na.omit(other_ds_dt)

agg_other_ds_dt <- aggregate(corr ~ symbol1+symbol2, data=other_ds_dt, FUN=aggFun)

count_other_ds_dt <- aggregate(corr ~ symbol1+symbol2, data=other_ds_dt, FUN=length)
colnames(count_other_ds_dt) <- c("symbol2_hk", "symbol1_hk", "count") # => reverse 2 and 1 for the other datasets


# oi_ds_dt$symbol1_hk <- paste0(oi_ds_dt$symbol1, ".x")
# oi_ds_dt$symbol2_hk <- paste0(oi_ds_dt$symbol2, ".y")
# agg_other_ds_dt$symbol1_hk <- paste0(agg_other_ds_dt$symbol2, ".x")
# agg_other_ds_dt$symbol2_hk <- paste0(agg_other_ds_dt$symbol1, ".y")

oi_ds_dt$symbol1_hk <- oi_ds_dt$symbol1
oi_ds_dt$symbol2_hk <- oi_ds_dt$symbol2
agg_other_ds_dt$symbol1_hk <- agg_other_ds_dt$symbol2 # => reverse 2 and 1 for the other datasets
agg_other_ds_dt$symbol2_hk <- agg_other_ds_dt$symbol1 


plot_dt <- rbind(oi_ds_dt[,c("symbol1_hk", "symbol2_hk", "corr")],
                 agg_other_ds_dt[,c("symbol1_hk", "symbol2_hk", "corr")])


plot_dt <- merge(plot_dt, count_other_ds_dt, all.x=TRUE, all.y=TRUE, by=c("symbol1_hk", "symbol2_hk"))

plot_dt$symbol1_hk <-factor(plot_dt$symbol1_hk, levels = symbols_level)
stopifnot(!is.na(plot_dt$symbol1_hk))
plot_dt$symbol2_hk <-factor(plot_dt$symbol2_hk, levels = symbols_level)
stopifnot(!is.na(plot_dt$symbol2_hk))

plot_dt$pair <- as.character(interaction(as.character(plot_dt$symbol1_hk), as.character(plot_dt$symbol2_hk)))
stopifnot(!duplicated(plot_dt$pair)) # ! if correct, no duplicated

plotTit <- paste0(hicds_oi, " - ", exprds_oi, " (upper left)")
subTit <- paste0(tad_oi, " (lower right: all other datasets)")
legTit <- "Expression correlation"
lowColor <- "blue"
highColor <- "red"

text_bl <- paste0("mean selected\n", round(mean(oi_ds_dt$corr), 2))
text_tr <- paste0("mean other\n", round(mean(agg_other_ds_dt$corr), 2))

mean_oi <- round(mean(oi_ds_dt$corr), 2)
mean_other <- round(mean(agg_other_ds_dt$corr), 2)

plot_dt$corr_rd <- round(plot_dt$corr, 2)
plot_dt$corr_rd_lab <- ifelse(is.na(plot_dt$count), plot_dt$corr_rd, paste0(plot_dt$corr_rd, "\n(",plot_dt$count, ")" ))

hm_p <- ggplot(data = plot_dt, aes(x=symbol1_hk, y=symbol2_hk, fill=corr)) + 
  ggtitle(paste0(plotTit), subtitle=paste0(subTit))+
  geom_tile()+
  geom_text(aes(label=corr_rd_lab),color = "black", size = 8, fontface="bold") +
  labs(fill=paste0(legTit))+
  scale_fill_gradient(low=lowColor, high=highColor) +
  theme(
    plot.title = element_text(size=18, hjust=0.5, face = "bold", family = "Hershey"),
    plot.subtitle = element_text(size=14, hjust=0.5, face = "italic", family = "Hershey"),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size=16, family = "Hershey"),
    axis.ticks = element_blank()
  ) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 2,
                               title.position = "top", title.hjust = 0.5))


hm_p_annot <- hm_p +  theme(
              legend.title = element_text(size=12, face="bold"),
              legend.position = c(1, 1),
              legend.justification = c(1,1), 
              legend.direction = "horizontal") + 
  # annotate("text", 
  #          x =1, y=1, 
  #          label = paste("atop(bold(mean~selected)==", mean_oi, ",\u03BC~bold(mean~other)==", mean_other, ")", sep = ""),
  #          # label = paste("\u2193"))
  #         parse=TRUE)

annotate("text", 
         x =1, y=c(1.2, 0.8), 
         # label = c( paste0("\u2191 mean selected = ", mean_oi ),
         #            paste0("\u2192 mean other = ", mean_other )), fontface="bold")
        label = c( paste0("mean selected = ", mean_oi, "\u2191" ),
                   paste0("mean other = ", mean_other, "\u2192" )), fontface="bold", size=5)


outFile <- file.path(outFolder, paste0(hicds_oi, "_", exprds_oi, "_", tad_oi, "_pairwiseCorr_cmpOtherDS.", plotType))
ggsave(plot = hm_p_annot, filename = outFile, height=myHeight, width=myWidth)





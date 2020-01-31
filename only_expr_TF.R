hicds <-  "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)

stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)


reg_file <- file.path("chea3_lung_TFs_processed.txt")
reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
cat(paste0("with Target Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)

reg_dt <- reg_dt[reg_dt$regSymbol %in% names(symb2entrez),]
cat(paste0("with Reg Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
reg_dt$regEntrezID <- symb2entrez[reg_dt$regSymbol]
reg_dt$regEntrezID <- as.character(reg_dt$regEntrezID)

all_reg_entrez <- unique(reg_dt$regEntrezID)

g2t_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2t_dt$entrezID <- as.character(g2t_dt$entrezID)

de_dt <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "1_runGeneDE", "DE_topTable.Rdata") ))

sum(all_reg_entrez %in% de_dt$genes)
length(all_reg_entrez)

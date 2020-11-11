setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

hicds <- "ENCSR489OCU_NCI-H460_40kb"

lg1_dt <- read.delim(paste0(hicds, "/genes2tad/all_genes_positions.txt"), col.names=c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
lg1_dt$entrezID <- as.character(lg1_dt$entrezID)
lg1_dt$symb <- entrez2symb[lg1_dt$entrezID]
lg1_t <- lg1_dt$region[lg1_dt$symb == "AKR1C1"]
lg1_dt[lg1_dt$region == lg1_t,]

lg1_tad_dt <- read.delim(paste0(hicds, "/genes2tad/all_assigned_regions.txt"), col.names=c( "chromo", "region", "start", "end"), stringsAsFactors = FALSE)
idx <- which(lg1_tad_dt$region == lg1_t)
lg1_tad_dt[(idx-1):(idx+1),]

# LG1:
#   chromo      region   start     end
# 15  chr10 chr10_TAD15 4600001 4880000
# 16  chr10 chr10_TAD16 4880001 5120000
# 17  chr10 chr10_TAD17 5120001 5480000

# LG2:
# chromo      region   start     end
# 14  chr10 chr10_TAD14 4520001 4920000
# 15  chr10 chr10_TAD15 4920001 5120000
# 16  chr10 chr10_TAD16 5120001 5200000

# ENCSR444WCZ_A549_40kb
# 20  chr10 chr10_TAD20 4840001 4920000
# 21  chr10 chr10_TAD21 4920001 5080000
# 22  chr10 chr10_TAD22 5080001 5240000

# ENCSR489OCU_NCI-H460_40kb
# 15  chr10 chr10_TAD15 4840001 4920000
# 16  chr10 chr10_TAD16 4920001 5120000
# 17  chr10 chr10_TAD17 5120001 5520000

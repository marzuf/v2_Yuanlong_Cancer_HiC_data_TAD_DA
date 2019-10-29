hicds="LG1_40kb"
exprds = "TCGAlusc_norm_lusc"

dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))


setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


mmp_genes <- entrez2symb[grepl("MMP", entrez2symb)]

dt[dt$entrezID %in% names(mmp_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(mmp_genes)]))

"chr10_TAD445" 
"chr11_TAD14" 
"chr11_TAD114"
"chr11_TAD389" 
"chr11_TAD390" 
"chr12_TAD198" 
"chr12_TAD493" 
"chr14_TAD13" 
"chr15_TAD10" 
"chr15_TAD34" 
"chr15_TAD44" 
"chr16_TAD16"  
"chr16_TAD152" 
"chr16_TAD163" 
"chr17_TAD109" 
"chr1_TAD6"   
"chr20_TAD121" 
"chr20_TAD156"
"chr22_TAD22" 
"chr7_TAD425" 
"chr8_TAD336" 


# Rscript plot_lolli_by_tad_dataset.R chr10_TAD445 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr11_TAD114 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr11_TAD14 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr11_TAD389 LG1_40kb TCGAlusc_norm_lusc 
# 
# Rscript plot_lolli_by_tad_dataset.R chr12_TAD198 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr11_TAD390 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr12_TAD493 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr14_TAD13 LG1_40kb TCGAlusc_norm_lusc 



# Rscript plot_lolli_by_tad_dataset.R chr15_TAD10 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr15_TAD34 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr15_TAD44 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr16_TAD16 LG1_40kb TCGAlusc_norm_lusc 
# 
# 
# Rscript plot_lolli_by_tad_dataset.R chr16_TAD152 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr16_TAD163 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr17_TAD109 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr1_TAD6 LG1_40kb TCGAlusc_norm_lusc 
# 
# 
# 
# Rscript plot_lolli_by_tad_dataset.R chr20_TAD121 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr20_TAD156 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr22_TAD22 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD425 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr8_TAD336 LG1_40kb TCGAlusc_norm_lusc 



six_genes <- entrez2symb[grepl("^SIX", entrez2symb)]

dt[dt$entrezID %in% names(six_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(six_genes)]))

"chr14_TAD143" 
"chr19_TAD154"
"chr2_TAD190"  
"chr2_TAD191"


# Rscript plot_lolli_by_tad_dataset.R chr14_TAD143 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr19_TAD154 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr2_TAD190 LG1_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr2_TAD191 LG1_40kb TCGAlusc_norm_lusc 


############################################################################################## KLK genes prostate

hicds="GSE118514_RWPE1_40kb"
exprds = "TCGAprad_norm_prad"

dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))


setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


klk_genes <- entrez2symb[grepl("KLK", entrez2symb)]

dt[dt$entrezID %in% names(klk_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(klk_genes)]))

unique(as.character(dt$region[dt$entrezID %in% names(klk_genes)]))
# [1] "chr19_TAD172" "chr19_TAD173" "chr4_TAD719" 

# Rscript plot_lolli_by_tad_dataset.R chr19_TAD173 GSE118514_RWPE1_40kb TCGAprad_norm_prad 
# Rscript plot_lolli_by_tad_dataset.R chr19_TAD172 GSE118514_RWPE1_40kb TCGAprad_norm_prad 
# Rscript plot_lolli_by_tad_dataset.R chr4_TAD719 GSE118514_RWPE1_40kb TCGAprad_norm_prad 

############################################################################################## IFIT genes liver
# Ch25H IFIT1 IFIT2 IFIT3 IFIT5 LIPA in liver cancer
# chr10_TAD316  

hicds="LI_40kb"
exprds = "TCGAlihc_norm_lihc"

dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))


setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


ifit_genes <- entrez2symb[grepl("IFIT", entrez2symb)]

dt[dt$entrezID %in% names(ifit_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(ifit_genes)]))

unique(as.character(dt$region[dt$entrezID %in% names(ifit_genes)]))
# [1] "chr10_TAD316" "chr11_TAD1"   "chr11_TAD4"   "chr11_TAD255" "chr13_TAD54" 
# [6] "chr6_TAD108"  "chr8_TAD209" 

# Rscript plot_lolli_by_tad_dataset.R chr10_TAD316 LI_40kb TCGAlihc_norm_lihc 
# Rscript plot_lolli_by_tad_dataset.R chr11_TAD1 LI_40kb TCGAlihc_norm_lihc 
# Rscript plot_lolli_by_tad_dataset.R chr11_TAD4 LI_40kb TCGAlihc_norm_lihc 
# Rscript plot_lolli_by_tad_dataset.R chr11_TAD255 LI_40kb TCGAlihc_norm_lihc 
# Rscript plot_lolli_by_tad_dataset.R chr13_TAD54 LI_40kb TCGAlihc_norm_lihc 
# Rscript plot_lolli_by_tad_dataset.R chr6_TAD108 LI_40kb TCGAlihc_norm_lihc 
# Rscript plot_lolli_by_tad_dataset.R chr8_TAD209 LI_40kb TCGAlihc_norm_lihc 


############################################################################################## GIMAP genes cancer
# ENCSR444WCZ_A549_40kb/ GIMAP 2 4 6  7
  # chr7_TAD553  
# "TCGAluad_norm_luad"

hicds="ENCSR444WCZ_A549_40kb"
exprds = "TCGAluad_norm_luad"
exprds = "TCGAlusc_norm_lusc"


dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))


setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


gimap_genes <- entrez2symb[grepl("GIMAP", entrez2symb)]

dt[dt$entrezID %in% names(gimap_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(gimap_genes)]))

unique(as.character(dt$region[dt$entrezID %in% names(gimap_genes)]))
# [1] "chr7_TAD552" "chr7_TAD553" "chr7_TAD554"


# Rscript plot_lolli_by_tad_dataset.R chr7_TAD552 ENCSR444WCZ_A549_40kb TCGAluad_norm_luad 
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD553 ENCSR444WCZ_A549_40kb TCGAluad_norm_luad 
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD554 ENCSR444WCZ_A549_40kb TCGAluad_norm_luad 

# Rscript plot_lolli_by_tad_dataset.R chr7_TAD552 ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD553 ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc 
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD554 ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc 


##########

data1 <- get(load("GO_SIGNIF_GENELEVEL_TADLEVEL_INTERSECTDIFF/tadPvalThresh0.01_genePvalThresh0.01/all_go_enrich_list.Rdata"))

data2 <- get(load("GO_SIGNIF_GENELEVEL_TADLEVEL_INTERSECTDIFF/tadPvalThresh0.05_genePvalThresh0.05/all_go_enrich_list.Rdata"))

length(data1[[1]][["tads_signif_genes"]]) < length(data2[[1]][["tads_signif_genes"]])
length(data1[[1]][["limma_signif_genes"]]) < length(data2[[1]][["limma_signif_genes"]])


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


############################################################################################## GIMAP genes cancer
# ENCSR444WCZ_A549_40kb/ GIMAP 2 4 6  7
# chr7_TAD553  
# "TCGAluad_norm_luad"

hicds="GSE99051_786_O_40kb"
exprds = "TCGAkich_norm_kich"


dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))


setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


gimap_genes <- entrez2symb[grepl("AKR1", entrez2symb)]

dt[dt$entrezID %in% names(gimap_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(gimap_genes)]))

unique(as.character(dt$region[dt$entrezID %in% names(gimap_genes)]))
# [1] "chr10_TAD23" "chr10_TAD24"  "chr10_TAD244" "chr11_TAD63"  "chr13_TAD107"
# [6] "chr14_TAD157" "chr18_TAD37"  "chr18_TAD243" "chr19_TAD124" "chr1_TAD171" 
# [11] "chr1_TAD536"  "chr1_TAD763"  "chr3_TAD279"  "chr7_TAD506"  "chr7_TAD519" 

# Rscript plot_lolli_by_tad_dataset.R chr10_TAD244 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_lolli_by_tad_dataset.R chr11_TAD63 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_lolli_by_tad_dataset.R chr13_TAD107 GSE99051_786_O_40kb TCGAkich_norm_kich

# Rscript plot_lolli_by_tad_dataset.R chr14_TAD157 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_lolli_by_tad_dataset.R chr18_TAD37 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_lolli_by_tad_dataset.R chr18_TAD243 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_lolli_by_tad_dataset.R chr19_TAD124 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_lolli_by_tad_dataset.R chr1_TAD171 GSE99051_786_O_40kb TCGAkich_norm_kich

# Rscript plot_lolli_by_tad_dataset.R chr1_TAD536 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_lolli_by_tad_dataset.R chr1_TAD763 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_lolli_by_tad_dataset.R chr3_TAD279 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD506 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD519 GSE99051_786_O_40kb TCGAkich_norm_kich




 

############################################################################################## C1Q

hicds="ENCSR312KHQ_SK-MEL-5_40kb"
exprds = "TCGAskcm_lowInf_highInf"


dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))


setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

fam_genes <- entrez2symb[grepl("^C1Q", entrez2symb)]

dt[dt$entrezID %in% names(fam_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(fam_genes)]))

# Rscript plot_lolli_by_tad_dataset.R chr10_TAD66 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr11_TAD193 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr11_TAD452 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr12_TAD189 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr13_TAD21 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr13_TAD24 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr16_TAD6 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr17_TAD24 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr1_TAD84 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr21_TAD35 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr22_TAD79 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr2_TAD437 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr17_TAD156 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr5_TAD140 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr4_TAD55 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr5_TAD613 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr17_TAD284 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf



############################################################################################## GIMAP

hicds="ENCSR312KHQ_SK-MEL-5_40kb"
exprds = "TCGAskcm_lowInf_highInf"

dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))

setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

fam_genes <- entrez2symb[grepl("^GIMAP", entrez2symb)]

dt[dt$entrezID %in% names(fam_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(fam_genes)]))

# Rscript plot_lolli_by_tad_dataset.R chr7_TAD569 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD568 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD567 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf



############################################################################################## S100

hicds="K562_40kb"
exprds = "TCGAlaml_wt_mutFLT3"

dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))

setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

fam_genes <- entrez2symb[grepl("^S100", entrez2symb)]

dt[dt$entrezID %in% names(fam_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(fam_genes)]))

# Rscript plot_lolli_by_tad_dataset.R chr11_TAD230 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_lolli_by_tad_dataset.R chr12_TAD204 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_lolli_by_tad_dataset.R chr1_TAD438 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_lolli_by_tad_dataset.R chr1_TAD118 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_lolli_by_tad_dataset.R chr1_TAD439 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_lolli_by_tad_dataset.R chr1_TAD444 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_lolli_by_tad_dataset.R chr1_TAD443 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD314 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD112 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_lolli_by_tad_dataset.R chr21_TAD112 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_lolli_by_tad_dataset.R chr5_TAD224 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_lolli_by_tad_dataset.R chr4_TAD27 K562_40kb TCGAlaml_wt_mutFLT3


############################################################################################## GIMAP

hicds="ENCSR862OGI_RPMI-7951_40kb"
exprds = "TCGAskcm_lowInf_highInf"

dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))

setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

fam_genes <- entrez2symb[grepl("^GIMAP", entrez2symb)]

dt[dt$entrezID %in% names(fam_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(fam_genes)]))

# Rscript plot_lolli_by_tad_dataset.R chr7_TAD537 ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD538 ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf
# Rscript plot_lolli_by_tad_dataset.R chr7_TAD539 ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf



############################################################################################## HOXC

hicds="ENCSR346DCU_LNCaP_40kb"
exprds = "TCGAprad_norm_prad"

dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))

setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

fam_genes <- entrez2symb[grepl("^HOXC", entrez2symb)]

dt[dt$entrezID %in% names(fam_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(fam_genes)]))

# Rscript plot_lolli_by_tad_dataset.R chr12_TAD220 ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad



 
############################################################################################## HOXC

hicds="GSE118514_22Rv1_40kb"
exprds = "TCGAprad_norm_prad"

dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))

setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

fam_genes <- entrez2symb[grepl("^HOXC", entrez2symb)]

dt[dt$entrezID %in% names(fam_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(fam_genes)]))

# Rscript plot_lolli_by_tad_dataset.R chr12_TAD195 GSE118514_22Rv1_40kb TCGAprad_norm_prad
# Rscript plot_lolli_by_tad_dataset.R chr12_TAD196 GSE118514_22Rv1_40kb TCGAprad_norm_prad



############################################################################################## AKR1C

hicds="LG2_40kb"
exprds = "TCGAluad_mutKRAS_mutEGFR"

dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))

setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

fam_genes <- entrez2symb[grepl("^AKR1C", entrez2symb)]

dt[dt$entrezID %in% names(fam_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(fam_genes)]))

# Rscript plot_lolli_by_tad_dataset.R chr10_TAD15 LG2_40kb TCGAluad_mutKRAS_mutEGFR
# Rscript plot_lolli_by_tad_dataset.R chr10_TAD16 LG2_40kb TCGAluad_mutKRAS_mutEGFR
# Rscript plot_lolli_by_tad_dataset.R chr10_TAD17 LG2_40kb TCGAluad_mutKRAS_mutEGFR


############################################################################################## SPRR

hicds="LG2_40kb"
exprds = "TCGAlusc_norm_lusc"

dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))

setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

fam_genes <- entrez2symb[grepl("^SPRR", entrez2symb)]

dt[dt$entrezID %in% names(fam_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(fam_genes)]))

# Rscript plot_lolli_by_tad_dataset.R chr1_TAD519 LG2_40kb TCGAlusc_norm_lusc
# Rscript plot_lolli_by_tad_dataset.R chr1_TAD520 LG2_40kb TCGAlusc_norm_lusc





############################################################################################## AKR1C

hicds="LG2_40kb"
exprds = "TCGAluad_nonsmoker_smoker"

dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))

setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

fam_genes <- entrez2symb[grepl("^AKR1C", entrez2symb)]

dt[dt$entrezID %in% names(fam_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(fam_genes)]))
# 
# Rscript plot_lolli_by_tad_dataset.R chr10_TAD15 LG2_40kb TCGAluad_nonsmoker_smoker
# Rscript plot_lolli_by_tad_dataset.R chr10_TAD16 LG2_40kb TCGAluad_nonsmoker_smoker
# Rscript plot_lolli_by_tad_dataset.R chr10_TAD17 LG2_40kb TCGAluad_nonsmoker_smoker




############################################################################################## HOXC

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

fam_genes <- entrez2symb[grepl("^HOXC", entrez2symb)]

dt[dt$entrezID %in% names(fam_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(fam_genes)]))

# Rscript plot_lolli_by_tad_dataset.R chr12_TAD194 GSE118514_RWPE1_40kb TCGAprad_norm_prad

 
############################################################################################## EPHB

hicds="ENCSR312KHQ_SK-MEL-5_40kb"
exprds = "TCGAskcm_lowInf_highInf"

dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))

setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

fam_genes <- entrez2symb[grepl("^EPHB", entrez2symb)]

dt[dt$entrezID %in% names(fam_genes),]

unique(as.character(dt$region[dt$entrezID %in% names(fam_genes)]))

Rscript plot_lolli_by_tad_dataset.R chr1_TAD84 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
Rscript plot_lolli_by_tad_dataset.R chr3_TAD521 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
Rscript plot_lolli_by_tad_dataset.R chr3_TAD722 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
Rscript plot_lolli_by_tad_dataset.R chr7_TAD375 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
Rscript plot_lolli_by_tad_dataset.R chr7_TAD538 ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf












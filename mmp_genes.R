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


  
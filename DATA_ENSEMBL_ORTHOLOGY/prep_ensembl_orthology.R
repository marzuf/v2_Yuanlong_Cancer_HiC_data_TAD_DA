 # Rscript prep_ensembl_orthology.R

id_file <- "ensembl_human_id_data.txt"
id_dt <- read.delim(id_file, header=TRUE, stringsAsFactors = FALSE)

orth_file <- "ensembl_orthology_data.txt"
all_orth_dt <- read.delim(orth_file, header=TRUE, stringsAsFactors = FALSE)
all_orth_dt$Gene.start..bp. <- NULL
all_orth_dt$Gene.end..bp. <- NULL
all_orth_dt$Gene.name <- NULL

all_na_ortho <- apply(all_orth_dt, 1, function(x) sum(is.na(x) | x == "") == 28 )  # 29 cols, 1 corresp. to human gene names/features
orth_dt <- all_orth_dt[!all_na_ortho,]


id_match_dt <- id_dt[,c("Gene.stable.ID", "EntrezGene.ID")]
id_match_dt <- unique(id_match_dt)

ensembl_orthology_prep_data <- merge(orth_dt, id_match_dt, all.x = TRUE, all.y = FALSE, by="Gene.stable.ID")

colnames(ensembl_orthology_prep_data)[colnames(ensembl_orthology_prep_data) == "EntrezGene.ID"] <- "entrezID"
ensembl_orthology_prep_data$entrezID <- as.character(ensembl_orthology_prep_data$entrezID)

outfile <- "ensembl_orthology_prep_data.Rdata"
save(ensembl_orthology_prep_data, file = outfile, version=2)
cat(paste0("... written: ", outfile, "\n"))




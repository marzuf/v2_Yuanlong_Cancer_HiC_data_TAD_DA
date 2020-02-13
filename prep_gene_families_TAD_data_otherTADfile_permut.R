library(dplyr)
library(reshape2)

# Rscript prep_gene_families_TAD_data_otherTADfile_permut.R ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad

startTime <- Sys.time()

# => copied "prep_gene_families_TAD_data.R" from /mnt/etemp/marie/TAD_DE_pipeline_v2_topRanking
# change the assigned_regions.txt for each of the dataset

SSHFS <- F
# SSHFS <- T
setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")
savePlot <- ifelse(SSHFS, FALSE, TRUE)

hgnc_geneFamilyFile <- file.path(setDir, "/mnt/ed4/marie/family_data_2/hgnc_entrez_family.txt")
ensembl_geneFamilyFile <- file.path(setDir, "/mnt/ed4/marie/family_data_2_ensembl/grch37_ensembl_family.csv")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args)==2)
curr_ds <- args[1]
hicds <- args[1]
exprds <- args[2]

#**************************************************************************************************
# PIPELINE SETTINGS
#**************************************************************************************************

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

cat(paste0("> Rscript prep_gene_families_TAD_data_otherTADfile_permut.R\n"))

# cat(paste0("setDir = ", setDir, "\n"))
# source("main_settings.R") # setDir is the main_settings not in run_settings
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

outFold <- file.path("PREP_GENE_FAMILIES_TAD_DATA_PERMUT", curr_ds, exprds)
dir.create(outFold, recursive = TRUE)

logFile <- file.path(outFold, paste0("prep_gene_families_logfile.txt"))
system(paste0("rm -f ", logFile))

caller <- "TopDom"

#**************************************************************************************************
# GENERAL DATA FROM THE PIPELINE
#**************************************************************************************************
nPermutRef=1
permut_dt <- get(load(paste0(hicds, "_", exprds, "_1000permut_permDT.Rdata" )))
gene2tadDT <- data.frame(entrezID = rownames(permut_dt), region = permut_dt[,nPermutRef], stringsAsFactors = FALSE)
# gene2tadDT_file <- file.path(curr_ds, "genes2tad", "all_genes_positions.txt")    
# gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

entrez2ensembl_file <- file.path(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_ENSEMBL/final_entrez2ensembl.txt")
entrez2ensembl_DT <- read.delim(entrez2ensembl_file, header=T, stringsAsFactors = F)
entrez2ensembl_DT$entrezID <- as.character(entrez2ensembl_DT$entrezID)

# file with coordinates of all regions  # => CHANGED HERE FOR THE DATASETS !!!

# TADpos_file <- paste0(setDir, "/mnt/ed4/marie/gene_data_final/consensus_TopDom_covThresh_r0.6_t80000_v0_w-1_final/genes2tad/all_assigned_regions.txt")    
# tadDT <- read.delim(TADpos_file, header=F, sep="\t", col.names = c("chromo", "region", "start", "end"), stringsAsFactors = F)
# TADpos_file <- file.path(curr_ds, "genes2tad", "all_assigned_regions.txt")    
# tadDT <- read.delim(TADpos_file, header=F, sep="\t", col.names = c("chromo", "region", "start", "end"), stringsAsFactors = F)


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TAKE ONLY TAD regions
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# tadDT <- tadDT[grepl("_TAD", tadDT$region),]
gene2tadDT <- gene2tadDT[grepl("_TAD", gene2tadDT$region),]
gene2tadDT$chromo <- gsub("(.+)_TAD.+", "\\1", gene2tadDT$region)
stopifnot(gene2tadDT$chromo %in% paste0("chr", 1:22))
entrez2ensembl_DT <- entrez2ensembl_DT[entrez2ensembl_DT$entrezID %in% gene2tadDT$entrezID,]

#**************************************************************************************************
# PREPARE THE HGNC FAMILY DATA
#**************************************************************************************************

hgnc_geneFamilyDT <- read.delim(hgnc_geneFamilyFile, col.names=c("entrezID", "family"), header = F, stringsAsFactors = F)
hgnc_geneFamilyDT$entrezID <- as.character(hgnc_geneFamilyDT$entrezID)
hgnc_geneFamilyDT$family_short <- unlist(sapply(hgnc_geneFamilyDT$family, function(x) strsplit(x, "\\|")[[1]][1] ))
# any(duplicated(hgnc_geneFamilyDT$entrezID))
# FALSE
stopifnot(!any(duplicated(hgnc_geneFamilyDT$entrezID)))
hgnc_geneFamilyDT <- hgnc_geneFamilyDT[hgnc_geneFamilyDT$family != "",]
hgnc_geneFamilyDT <- na.omit(hgnc_geneFamilyDT)
length(hgnc_geneFamilyDT$entrezID)
# 19224
length(unique(hgnc_geneFamilyDT$entrezID))
# 19224

hgnc_geneFamilyDT <- hgnc_geneFamilyDT[,c("entrezID", "family", "family_short")]
colnames(hgnc_geneFamilyDT) <- c("entrezID", "hgnc_family", "hgnc_family_short")

gene2tadDT_hgncFamily <- left_join(gene2tadDT, hgnc_geneFamilyDT, by="entrezID")
# gene2tadDT_hgncFamily$hgnc_family[is.na(gene2tadDT_hgncFamily$hgnc_family)] <- paste0(gene2tadDT_hgncFamily$region[is.na(gene2tadDT_hgncFamily$hgnc_family)], "_", 
#                                                                                       gene2tadDT_hgncFamily$entrezID[is.na(gene2tadDT_hgncFamily$hgnc_family)])
# 
# gene2tadDT_hgncFamily$hgnc_family_short[is.na(gene2tadDT_hgncFamily$hgnc_family_short)] <- paste0(gene2tadDT_hgncFamily$region[is.na(gene2tadDT_hgncFamily$hgnc_family_short)], "_", 
#                                                                                                   gene2tadDT_hgncFamily$entrezID[is.na(gene2tadDT_hgncFamily$hgnc_family_short)])

#**************************************************************************************************
# PREPARE THE ENSEMBL FAMILY DATA
#**************************************************************************************************
ensembl_geneFamilyDT <- read.delim(ensembl_geneFamilyFile, sep=",", header = T, stringsAsFactors = F)
ensembl_geneFamilyDT$Transcript.stable.ID <- NULL
ensembl_geneFamilyDT <- unique(ensembl_geneFamilyDT)
stopifnot(nrow(ensembl_geneFamilyDT) > 0)
ensembl_geneFamilyDT$Chromosome.scaffold.name <- paste0("chr", ensembl_geneFamilyDT$Chromosome.scaffold.name)
ensembl_geneFamilyDT <- ensembl_geneFamilyDT[ensembl_geneFamilyDT$Chromosome.scaffold.name %in% gene2tadDT$chromo,]
ensembl_geneFamilyDT <- unique(ensembl_geneFamilyDT)
stopifnot(nrow(ensembl_geneFamilyDT) > 0)
ensembl_geneFamilyDT <- ensembl_geneFamilyDT[ensembl_geneFamilyDT$Ensembl.Protein.Family.ID.s. != "",]
ensembl_geneFamilyDT <- unique(ensembl_geneFamilyDT)
ensembl_geneFamilyDT$family_short <- gsub("(.+)_.+", "\\1", ensembl_geneFamilyDT$Ensembl.Protein.Family.ID.s.)
# RETRIEVE THE entrezID FROM THE ensemblID
ensembl_geneFamilyDT <- left_join(ensembl_geneFamilyDT, entrez2ensembl_DT, by=c("Gene.stable.ID" = "ensemblID"))
colnames(ensembl_geneFamilyDT)[colnames(ensembl_geneFamilyDT) == "Ensembl.Protein.Family.ID.s."] <- "family"
colnames(ensembl_geneFamilyDT)[colnames(ensembl_geneFamilyDT) == "Ensembl.Family.Description"] <- "family_desc"

ensembl_geneFamilyDT <- ensembl_geneFamilyDT[,c("entrezID", "family", "family_short", "family_desc")]
ensembl_geneFamilyDT <- na.omit(ensembl_geneFamilyDT)
ensembl_geneFamilyDT <- unique(ensembl_geneFamilyDT)
any(duplicated(ensembl_geneFamilyDT$entrezID))
ensembl_geneFamilyDT[ensembl_geneFamilyDT$entrezID == ensembl_geneFamilyDT$entrezID[which(duplicated(ensembl_geneFamilyDT$entrezID))[1]],]
length(ensembl_geneFamilyDT$entrezID)
# 20402
length(unique(ensembl_geneFamilyDT$entrezID))
# 19372

ensembl_geneFamilyDT <- ensembl_geneFamilyDT[,c("entrezID", "family", "family_short")]
colnames(ensembl_geneFamilyDT) <- c("entrezID", "ensembl_family", "ensembl_family_short")

gene2tadDT_ensemblFamily <- left_join(gene2tadDT, ensembl_geneFamilyDT, by="entrezID")
# gene2tadDT_ensemblFamily$ensembl_family[is.na(gene2tadDT_ensemblFamily$ensembl_family)] <- paste0(gene2tadDT_ensemblFamily$region[is.na(gene2tadDT_ensemblFamily$ensembl_family)], "_", 
#                                                                                                   gene2tadDT_ensemblFamily$entrezID[is.na(gene2tadDT_ensemblFamily$ensembl_family)])
# 
# gene2tadDT_ensemblFamily$ensembl_family_short[is.na(gene2tadDT_ensemblFamily$ensembl_family_short)] <- paste0(gene2tadDT_ensemblFamily$region[is.na(gene2tadDT_ensemblFamily$ensembl_family_short)], "_", 
#                                                                                                               gene2tadDT_ensemblFamily$entrezID[is.na(gene2tadDT_ensemblFamily$ensembl_family_short)])

#**************************************************************************************************
#**************************************************************************************************
#**************************************************************************************************

head(gene2tadDT_hgncFamily)
hgncDT <- gene2tadDT_hgncFamily
txt <- paste0("... number of genes: ", nrow(hgncDT), "\n")
printAndLog(txt, logFile)
# cat("... number of genes with HGNC families: ", sum(!grepl("_TAD", hgncDT$hgnc_family)), "\n")
# cat("... number of genes with HGNC super-families: ", sum(!grepl("_TAD", hgncDT$hgnc_family_short)), "\n")
# hgncDT <- hgncDT[!grepl("_TAD", hgncDT$hgnc_family),]
txt <- paste0("... number of genes with HGNC families: ", sum(!is.na(hgncDT$hgnc_family)), "\n")
printAndLog(txt, logFile)
txt <- paste0("... number of genes with HGNC super-families: ", sum(!is.na(hgncDT$hgnc_family_short)), "\n")
printAndLog(txt, logFile)

hgncDT <- hgncDT[!is.na(hgncDT$hgnc_family),]

txt <- paste0("... number unique HGNC families: ", length(unique(hgncDT$hgnc_family)), "\n")
printAndLog(txt, logFile)
txt <- paste0("... number unique HGNC super-families: ", length(unique(hgncDT$hgnc_family_short)), "\n")
printAndLog(txt, logFile)
txt <- paste0("... summary gene HGNC families:\n")
printAndLog(txt, logFile)
sink(file=logFile, append=T)
summary(as.numeric(table(hgncDT$hgnc_family)))
sink()
txt <- paste0("\n")
printAndLog(txt, logFile)
txt <- paste0("... summary gene HGNC super-families:\n")
printAndLog(txt, logFile)
sink(file=logFile, append=T)
summary(as.numeric(table(hgncDT$hgnc_family_short)))
sink()
txt <- paste0("\n")
printAndLog(txt, logFile)

head(gene2tadDT_ensemblFamily)
ensemblDT <- gene2tadDT_ensemblFamily
txt <- paste0("... number of genes: ", nrow(ensemblDT), "\n")
printAndLog(txt, logFile)
# cat("... number of genes with ensembl families: ", sum(!grepl("_TAD", ensemblDT$ensembl_family)), "\n")
# cat("... number of genes with ensembl super-families: ", sum(!grepl("_TAD", ensemblDT$ensembl_family_short)), "\n")
# ensemblDT <- ensemblDT[!grepl("_TAD", ensemblDT$ensembl_family),]
txt <- paste0("... number of genes with ensembl families: ", sum(!is.na(ensemblDT$ensembl_family)), "\n")
printAndLog(txt, logFile)
txt <- paste0("... number of genes with ensembl super-families: ", sum(!is.na(ensemblDT$ensembl_family_short)), "\n")
printAndLog(txt, logFile)
ensemblDT <- ensemblDT[!is.na(ensemblDT$ensembl_family),]

txt <- paste0("... number unique ensembl families: ", length(unique(ensemblDT$ensembl_family)), "\n")
printAndLog(txt, logFile)
txt <- paste0("... number unique ensembl super-families: ", length(unique(ensemblDT$ensembl_family_short)), "\n")
printAndLog(txt, logFile)
txt <- paste0("... summary gene ensembl families:\n")
printAndLog(txt, logFile)
sink(file=logFile, append=T)
summary(as.numeric(table(ensemblDT$ensembl_family)))
sink()
txt <- paste0("\n")
printAndLog(txt, logFile)
txt <- paste0("... summary gene ensembl super-families:\n")
printAndLog(txt, logFile)
sink(file=logFile, append=T)
summary(as.numeric(table(ensemblDT$ensembl_family_short)))
sink()
txt <- paste0("\n")
printAndLog(txt, logFile)

#**************************************************************************************************
#************************************************************************************************** SAVE THE DATA
#**************************************************************************************************

hgnc_entrezID_family_TAD_DT <- hgncDT
outFile <- file.path(outFold, "hgnc_entrezID_family_TAD_DT.Rdata")
save(hgnc_entrezID_family_TAD_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

ensembl_entrezID_family_TAD_DT <- ensemblDT
outFile <- file.path(outFold, "ensembl_entrezID_family_TAD_DT.Rdata")
save(ensembl_entrezID_family_TAD_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))









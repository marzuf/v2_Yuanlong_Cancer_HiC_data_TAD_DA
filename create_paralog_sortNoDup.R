startTime <- Sys.time()

cat(paste0("> Rscript create_paralog_sortNoDup.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)



# Rscript create_paralog_sortNoDup.R 

# 1 file for all datasets !

### UPDATE sortNoDup 30.06.2018
# -> sort rows of the gene coord. table to ensure alphabetical order of the genes !!
#    so that after melt the 1st gene will always be the 1st in alphabetical order


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))


outFold <- file.path("CREATE_PARALOG_SORTNODUP")
dir.create(outFold, recursive=TRUE)

################################################################################################################################################
######################################################### PREP PARALOG DATA
################################################################################################################################################

id_dt <- read.delim("ensembl_entrez_biomart.txt", sep=",", stringsAsFactors=FALSE, header=TRUE)

paralog_dt <- read.delim("ensembl_paralogs_biomart.txt", sep=",", stringsAsFactors=FALSE, header=TRUE)
paralog_dt <- paralog_dt[!paralog_dt$Human.paralogue.gene.stable.ID=="",]

par_dt_sub <- paralog_dt[, c("Gene.stable.ID", "Human.paralogue.gene.stable.ID" )]
id_dt_sub <- id_dt[,c("Gene.stable.ID", "EntrezGene.ID")]
id_dt_sub$EntrezGene.ID <- as.character(id_dt_sub$EntrezGene.ID)

out_dt_tmp <- merge(par_dt_sub, id_dt_sub, by="Gene.stable.ID")

colnames(id_dt_sub)[colnames(id_dt_sub) == "Gene.stable.ID"] <- "Human.paralogue.gene.stable.ID"
colnames(id_dt_sub)[colnames(id_dt_sub) == "EntrezGene.ID"] <- "Human.paralogue.entrezGene.ID"

out_dt <- merge(out_dt_tmp, id_dt_sub, by="Human.paralogue.gene.stable.ID")
entrez_paraDT <- out_dt[,c("EntrezGene.ID", "Human.paralogue.entrezGene.ID")]
entrez_paraDT <- na.omit(entrez_paraDT)
nrow(entrez_paraDT)
entrez_paraDT$paralogs <- "paralogs"

entrez_paraDT$gene1 <- as.character(pmin(as.character(entrez_paraDT$EntrezGene), as.character(entrez_paraDT$Human.paralogue.entrezGene.ID)))
entrez_paraDT$gene2 <- as.character(pmax(as.character(entrez_paraDT$EntrezGene), as.character(entrez_paraDT$Human.paralogue.entrezGene.ID)))

 # ENSG00000005075       102113565     102119354     -1          5439    POLR2J
 # ENSG00000168255       102178365     102213103     -1          5439   POLR2J3

# ENSG00000005075                ENSG00000267645          
# ENSG00000005075                ENSG00000168255          
# ENSG00000005075                ENSG00000228049          
# ENSG00000005075                ENSG00000270249          
# ENSG00000168255                ENSG00000267645          
# ENSG00000168255                ENSG00000228049          
# ENSG00000168255                ENSG00000005075          
# ENSG00000168255                ENSG00000270249          

# => this will lead to gene1 == gene2
entrez_paraDT <- entrez_paraDT[entrez_paraDT$gene2 != entrez_paraDT$gene1,]
stopifnot(entrez_paraDT$gene1 < entrez_paraDT$gene2)

all_paralog_pairs <- entrez_paraDT[,c("gene1", "gene2", "paralogs")]
all_paralog_pairs <- unique(all_paralog_pairs)

outFile <- file.path(outFold, paste0("all_paralog_pairs.Rdata"))
save(all_paralog_pairs, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


  
######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

startTime <- Sys.time()
cat(paste0("> Rscript create_sameTADpermutG2t_sortNoDup_otherTADfile.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

# Rscript create_sameTADpermutG2t_sortNoDup_otherTADfile.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad


### UPDATE sortNoDup 30.06.2018
# -> sort rows of the gene coord. table to ensure alphabetical order of the genes !!
#    so that after melt the 1st gene will always be the 1st in alphabetical order
# -> add as.character() in apply !!! use "gene1" etc. instead of 1 index in apply 


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2) 

ds <- args[1]
hicds <- args[1]
exprds <-args[2]

outFold <- file.path("CREATE_SAME_TAD_PERMUTG2T_SORTNODUP_PERM", ds, exprds)
# gene2tadDT_file <- file.path(ds, "genes2tad", "all_genes_positions.txt")

dir.create(outFold, recursive=TRUE)
# stopifnot(file.exists(gene2tadDT_file))

nPermutRef <- 1

permut_dt <- get(load(paste0(hicds, "_", exprds, "_1000permut_permDT.Rdata" )))
gene2tadDT <- data.frame(entrezID = rownames(permut_dt), region = permut_dt[,nPermutRef], stringsAsFactors = FALSE)
# gene2tadDT <- read.delim(gene2tadDT_file, header = F, stringsAsFactors = F, col.names = c("entrezID", "chromo", "start", "end", "region"))
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
stopifnot(!any(duplicated(gene2tadDT$entrezID)))
gene2tadDT <- gene2tadDT[grep("_TAD", gene2tadDT$region),]

all_tads <- unique(gene2tadDT$region)

all_TAD_pairs <- foreach(tad = all_tads, .combine='rbind') %dopar% {
  tad_g2t_dt <- gene2tadDT[gene2tadDT$region == tad,]
  if(nrow(tad_g2t_dt) == 1) return(NULL)
  # UPDATE 30.06.2018 -> ENSURE AS.CHARACTER + ALPHABETICAL ORDER !!!
  tad_g2t_dt$entrezID <- as.character(tad_g2t_dt$entrezID)
  tad_g2t_dt <- tad_g2t_dt[order(tad_g2t_dt$entrezID),]
  tadDT <- as.data.frame(t(combn(tad_g2t_dt$entrezID, m=2)))
  colnames(tadDT) <- c("gene1", "gene2")
  tadDT$region <- tad
  tadDT$gene1 <- as.character(tadDT$gene1)
  tadDT$gene2 <- as.character(tadDT$gene2)
  stopifnot(tadDT$gene1 < tadDT$gene2)
  tadDT
}

outFile <- file.path(outFold, "all_TAD_pairs.Rdata")
save(all_TAD_pairs, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

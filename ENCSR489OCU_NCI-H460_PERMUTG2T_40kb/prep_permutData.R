# RANDOM DATA: TAKE THE G2T ASSIGNMENT OF THE PERMUT

# prepare the new all_genes_positions.txt
# file in the genes2tad
# the all_assigned_regions.txt is retrieved and kept unchanged from true data

# Rscript prep_permutData.R

require(doMC)
require(foreach)
registerDoMC(40)

runFolder <- ".."
hicds <- "ENCSR489OCU_NCI-H460_40kb"

outFolder <- "genes2tad"
dir.create(outFolder, recursive=T)

permRef <- 1

perm_dt <- get(load("../ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_1000permut_permDT.Rdata"))

# NB to be faster (saved from previous script)
#> perm_dt_1 <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/5_runPermutationsMedian/permutationsDT.Rdata"))
#> perm_dt_2 <- get(load("ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/5_runPermutationsMedian/permutationsDT.Rdata"))                                                                                                                                                              
#> perm_dt_2 <- get(load("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_1000permut_permDT.Rdata"))
#> stopifnot(perm_dt_1[,1] ==perm_dt_2[,1])

# copy the assigned regions
cmd <- paste("cp", file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt"), file.path(outFolder, "all_assigned_regions.txt" ))
cat(paste0("> ", cmd, "\n"))
system(cmd)

# 
g2t_file <- file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo","start", "end", "region"))
g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
g2t_dt$size <- g2t_dt$end-g2t_dt$start+1
stopifnot(!duplicated(g2t_dt$entrezID))
geneSize <- setNames(g2t_dt$size, g2t_dt$entrezID)
rm(g2t_dt)
#
tad_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt")
tad_dt <- read.delim(tad_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo","region", "start", "end"))


# new gene2tad assignment
g2t_permut <- data.frame(
  entrezID = as.character(rownames(perm_dt)),
  region = perm_dt[,permRef],
  stringsAsFactors = FALSE
)



# it will be a bit weird, because for the moment I don't change
# the start and end positions of the gene, only the chromo
# so a TAD will contain genes that have non adjacent positions...
# the positions of the genes are needed for the 5sameNbr_runPermutationsCorr.R
# first idea: (I infer chromo from the new random region, but keep start and end positions)
# new idea: sample start position within the TAD
# and keep true gene size

stopifnot(g2t_permut$region  %in% tad_dt$region)

rd_g2t_tad_dt <- merge(g2t_permut, tad_dt, by="region")
colnames(rd_g2t_tad_dt)[colnames(rd_g2t_tad_dt) == "start"] <- "region_start"
colnames(rd_g2t_tad_dt)[colnames(rd_g2t_tad_dt) == "end"] <- "region_end"

rd_g2t_tad_dt$rd_gene_start <- foreach(i=1:nrow(rd_g2t_tad_dt), .combine='c') %dopar% {
  sample(rd_g2t_tad_dt$region_start[i]:rd_g2t_tad_dt$region_end[i],1)
}

rd_g2t_tad_dt$rd_gene_end <- foreach(i=1:nrow(rd_g2t_tad_dt), .combine='c') %dopar% {
  stopifnot(rd_g2t_tad_dt$entrezID[i] %in% names(geneSize))
  rd_g2t_tad_dt$rd_gene_start[i] + geneSize[paste0(rd_g2t_tad_dt$entrezID[i])]
}

stopifnot(rd_g2t_tad_dt$rd_gene_end > rd_g2t_tad_dt$rd_gene_start)

stopifnot(setequal(rd_g2t_tad_dt$entrezID, g2t_permut$entrezID))

outFile <- file.path(outFolder, "all_genes_positions.txt")

cat("OUTFILE:", outFile, "\n")

rd_g2t_tad_dt <- rd_g2t_tad_dt[order(rd_g2t_tad_dt$chromo, rd_g2t_tad_dt$rd_gene_start, rd_g2t_tad_dt$rd_gene_end),]

write.table(rd_g2t_tad_dt[,c("entrezID", "chromo", "rd_gene_start", "rd_gene_end", "region" )], col.names = F, row.names = F, quote=F, append=F, sep="\t", file = outFile)
cat(paste0("... written: ", outFile, "\n"))







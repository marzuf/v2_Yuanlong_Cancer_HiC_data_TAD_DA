setDir <- "/media/electron"
setDir <- ""

# Rscript revision_expressionData.R
script_name="revision_expressionData.R"
startTime <- Sys.time()
cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
registerDoMC(40)

minGenes <- 3

outFolder <- "REVISION_EXPRESSIONDATA"
dir.create(outFolder, recursive = TRUE)

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


expr_var <- "FPKM"
symb_aggFun <- "mean"  # aggregate expression data to have unique gene symbol
rep_aggFun <- "mean"  # aggregate expression across replicates

# all_mrna_data=get(load("sub_expression_allDS.Rdata"))

all_mrna_data <-get(load(file.path(setDir,
                                   "/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/1.Data/exp_CELL_LINEs/exp_info.Rdata")))

all_ds <-  c(
  "LI_40kb"="Compendium_LI",
  "GSE105381_HepG2_40kb"="mega_GSE105381_HepG2",
  "LG1_40kb" ="Compendium_LG1",
  "ENCSR444WCZ_A549_40kb"="mega_ENCSR444WCZ_A549",
  "LG2_40kb"="Compendium_LG2",
  "ENCSR489OCU_NCI-H460_40kb"="mega_ENCSR489OCU_NCI-H460",
  "GSE118514_RWPE1_40kb"="mega_GSE118514_RWPE1",
  "ENCSR346DCU_LNCaP_40kb"="mega_ENCSR346DCU_LNCaP",
  "GSE118514_22Rv1_40kb"="GSE118514_22Rv1"
)

all_pairs <- c(
  file.path("LI_40kb","GSE105381_HepG2_40kb", "TCGAlihc_norm_lihc"),
  file.path("LG1_40kb" ,"ENCSR444WCZ_A549_40kb", "TCGAluad_norm_luad"),
  file.path("LG2_40kb" ,"ENCSR444WCZ_A549_40kb" ,"TCGAluad_norm_luad"),
  file.path("LG1_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("LG2_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("GSE118514_RWPE1_40kb", "ENCSR346DCU_LNCaP_40kb", "TCGAprad_norm_prad"),
  file.path("GSE118514_RWPE1_40kb", "GSE118514_22Rv1_40kb", "TCGAprad_norm_prad")
)


final_dt_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))
final_dt$regionID <- file.path(final_dt$hicds, final_dt$exprds, final_dt$region)
final_dt$signif_lab <- ifelse(final_dt$adjPvalComb <= 0.01, "signif.", "not signif.")

final_dt <- final_dt[final_dt$hicds %in% names(all_ds) & final_dt$exprds %in% basename(all_pairs),]
stopifnot(nrow(final_dt) > 0)
final_dt$datasets <- file.path(final_dt$hicds, final_dt$exprds)

final_dt$mRNA_mean <- NA
final_dt$mRNA_nGenes <- NA
final_dt$tad_nGenes <- NA

ds=unique(final_dt$datasets)[1]
ds="ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad"

final_mRNA_dt <- foreach(ds = unique(final_dt$datasets), .combine='rbind') %do% {
  
  cat(paste0("> START : ", ds, "\n"))
  
  exprds <- basename(ds)
  hicds <- dirname(ds)
  stopifnot(hicds %in% names(all_ds))
  
  ds <- as.character(all_ds[hicds])
  
  stopifnot(ds %in% names(all_mrna_data))
  
  ds_mrna_data <- all_mrna_data[[ds]]
  
  sub_final_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds,]
  stopifnot(nrow(sub_final_dt) > 0)
  
  expr_values <- lapply(ds_mrna_data, function(x) {
    stopifnot("external_gene_name" %in% colnames(x))
    stopifnot(paste0(expr_var) %in% colnames(x))
    # stopifnot(!duplicated(x[["external_gene_name"]])) ### NOT TRUE
    tmp_dt <- aggregate(as.formula(paste0(expr_var, "~external_gene_name")), data=x, FUN=symb_aggFun, na.rm=TRUE)
    stopifnot(!duplicated(tmp_dt[["external_gene_name"]])) ### NOT TRUE
    setNames(tmp_dt[[paste0(expr_var)]], tmp_dt[["external_gene_name"]])
  })
  # first check that they have the same gene list
  stopifnot(length(unique(lapply(expr_values, names))) == 1)
  # as they are the same, I can sum
  expr_values_agg <- apply( do.call(rbind, expr_values), 2,rep_aggFun)
  stopifnot(names(expr_values_agg) == names(expr_values[[1]]))
  
  # some checks
  tmp <- sort(expr_values_agg, decreasing = TRUE)
  stopifnot(tmp[1] == do.call(rep_aggFun, list(unlist(lapply(expr_values, function(x)x[[names(tmp)[1]]])))))
  stopifnot(tmp[3] == do.call(rep_aggFun, list(unlist(lapply(expr_values, function(x)x[[names(tmp)[3]]])))))  
  

  i=1
  
  stopifnot(nrow(sub_final_dt) >= 1)
  
  all_tads_dt <- foreach(i = 1:nrow(sub_final_dt), .combine='rbind') %dopar% {
    
    region_symbols <- as.character(unlist(strsplit(sub_final_dt$region_genes[i], split=",")))
    stopifnot(!duplicated(region_symbols))
    stopifnot(length(region_symbols) >= 3)
    
    #expr_values_agg[paste0(region_symbols)] is less clean because return NA if not any
    tad_mrna <- expr_values_agg[names(expr_values_agg) %in% paste0(region_symbols)] 
    nGenes <- sum(region_symbols %in% names(expr_values_agg))
    stopifnot(nGenes == length(tad_mrna))
    
    sub_final_dt$mRNA_mean[i] <- mean(tad_mrna)
    sub_final_dt$mRNA_nGenes[i] <- nGenes
    sub_final_dt$tad_nGenes[i] <- length(region_symbols) 
    
    sub_final_dt[i,]
  }
  outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_all_tads_dt.Rdata"))
  save(all_tads_dt, file=outFile, version=2)
  cat(paste0("written: ", outFile, "\n"))
  
  all_tads_dt
}

outFile <- file.path(outFolder,"final_mRNA_dt.Rdata")
save(final_mRNA_dt, file=outFile, version=2)
cat(paste0("written: ", outFile, "\n"))

######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

# tmp <- list(c("A"=1,"B"=2,"C"=3), c("A"=3,"B"=6,"C"=9), c("A"=2, "B"=1, "C" =3))
# Reduce(rep_aggFun, tmp)  ### !!! not working
# Reduce(`+`, tmp)/length(tmp)
# colMeans(do.call(rbind, tmp))
# all(colMeans(do.call(rbind, tmp)) == Reduce(`+`, tmp)/length(tmp))
# 
# apply( do.call(rbind, tmp), 2,rep_aggFun)
# all(apply( do.call(rbind, tmp), 2,rep_aggFun) == Reduce(`+`, tmp)/length(tmp))

# work only if length >= 2
# all(sapply(lapply(expr_values[2:length(expr_values)], names), FUN = identical, names(expr_values[[1]])))
# all(sapply(lapply(tmp[2:length(tmp)], names), FUN = identical, names(tmp[[1]])))
# > tmp <- list(c("A"=1,"B"=2,"C"=3), c("A"=3,"B"=6,"C"=9))
# > unique(lapply(tmp, names))
# [[1]]
# [1] "A" "B" "C"
# 
# > Reduce(`+`, tmp)
# A  B  C 
# 4  8 12 

# > tmp <- list(c(1,2,3), c(1,2,3))
# > unique(tmp)
# [[1]]
# [1] 1 2 3
# 
# > tmp <- list(c(1,2,3), c(1,2,4))
# > unique(tmp)
# [[1]]
# [1] 1 2 3
# 
# [[2]]
# [1] 1 2 4
# 
# > 
  
# > sub_data2=x[["mega_GSE105318_DLD1"]]
# > str(sub_data2)
# List of 2
# $ rep1:Classes ‘data.table’ and 'data.frame':  54849 obs. of  9 variables:
#   ..$ ensembl_gene_id   : chr [1:54849] "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" "ENSG00000000938" ...
# ..$ TPM               : num [1:54849] 79.87 2.28 17.89 0 0.17 ...
# ..$ FPKM              : num [1:54849] 70.91 2.03 15.88 0 0.15 ...
# ..$ expected_count    : num [1:54849] 499 81.7 256.3 0 2 ...
# ..$ external_gene_name: chr [1:54849] "DPM1" "SCYL3" "C1orf112" "FGR" ...
# ..$ chromosome_name   : num [1:54849] 20 1 1 1 1 6 6 6 1 1 ...
# ..$ start_position    : int [1:54849] 49551404 169818772 169631245 27938575 196621008 143815948 53362139 41040684 24683489 24742284 ...
# ..$ end_position      : int [1:54849] 49575092 169863408 169823221 27961788 196716634 143832827 53481768 41067715 24743424 24799466 ...
# ..$ gene_biotype      : chr [1:54849] "protein_coding" "protein_coding" "protein_coding" "protein_coding" ...
# ..- attr(*, ".internal.selfref")=<externalptr> 
#   ..- attr(*, "sorted")= chr "ensembl_gene_id"
# $ rep2:Classes ‘data.table’ and 'data.frame':  54849 obs. of  9 variables:
#   ..$ ensembl_gene_id   : chr [1:54849] "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" "ENSG00000000938" ...
# ..$ TPM               : num [1:54849] 88.64 2.61 16.89 0 0.2 ...
# ..$ FPKM              : num [1:54849] 75.12 2.21 14.32 0 0.17 ...
# ..$ expected_count    : num [1:54849] 627 105 279 0 2 ...
# ..$ external_gene_name: chr [1:54849] "DPM1" "SCYL3" "C1orf112" "FGR" ...
# ..$ chromosome_name   : num [1:54849] 20 1 1 1 1 6 6 6 1 1 ...
# ..$ start_position    : int [1:54849] 49551404 169818772 169631245 27938575 196621008 143815948 53362139 41040684 24683489 24742284 ...
# ..$ end_position      : int [1:54849] 49575092 169863408 169823221 27961788 196716634 143832827 53481768 41067715 24743424 24799466 ...
# ..$ gene_biotype      : chr [1:54849] "protein_coding" "protein_coding" "protein_coding" "protein_coding" ...
# ..- attr(*, ".internal.selfref")=<externalptr> 
#   
  
  # -> for each TAD
  # -> how many genes available in the  rna expr
  # -> mean FPKM
  # -> then I can do the same as for norm and tumor
# Danielle prefer TPM, but I use FPKM
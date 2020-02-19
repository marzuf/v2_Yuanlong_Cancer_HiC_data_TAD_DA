options(scipen=100)

SSHFS=F

setDir <- "/media/electron"
setDir <- ""

# Rscript cmp_permut_geneRank_tadRank.R


script_name <- "cmp_permut_geneRank_tadRank.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(RColorBrewer)
require(reshape2)
require(ggplot2)
# require(gplots)
registerDoMC(4)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


plotType <- "png"
myWidth <- myHeight <- 400

outFolder <- "CMP_PERMUT_GENERANK_TADRANK"
dir.create(outFolder, recursive = TRUE)

obs_rank_dt <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
obs_rank_dt$dataset <- file.path(obs_rank_dt$hicds, obs_rank_dt$exprds)

permut_rank_dt <- get(load("GENE_RANK_TAD_RANK_PERMUT/all_gene_tad_signif_dt.Rdata"))
permut_rank_dt$dataset <- file.path(permut_rank_dt$hicds, permut_rank_dt$exprds)
permut_rank_dt$hicds_lab <- gsub("(.+)_.+_40kb", "\\1_40kb", permut_rank_dt$hicds)
permut_rank_dt$dataset_lab <- file.path(permut_rank_dt$hicds_lab, permut_rank_dt$exprds)


ds_withPerm <- intersect(permut_rank_dt$dataset_lab, obs_rank_dt$dataset)
plot_ds = ds_withPerm[1]

for(plot_ds in ds_withPerm) {
  
  sub_obs_dt <- obs_rank_dt[obs_rank_dt$dataset == plot_ds,]
  stopifnot(nrow(sub_obs_dt)>0)
  
  sub_perm_dt <- permut_rank_dt[permut_rank_dt$dataset_lab == plot_ds,]
  all_perms <- unique(sub_perm_dt$dataset)
  rd=all_perms[1]
  
  for(rd in all_perms) {
    
    rd_type <- gsub(".+_(.+?_40kb)", "\\1", dirname(rd))
    
    rd_sub_perm_dt <- sub_perm_dt[sub_perm_dt$dataset == rd,]
    stopifnot(nrow(rd_sub_perm_dt) > 0)
    
    stopifnot(!duplicated(rd_sub_perm_dt$entrezID))
    stopifnot(!duplicated(sub_obs_dt$entrezID))
    
    merge_dt <- merge(sub_obs_dt[,c("entrezID", "tad_rank")],
                      rd_sub_perm_dt[,c("entrezID", "tad_rank")], 
                      by="entrezID",
                      all.x=FALSE, all.y=FALSE,
                      suffixes = c("_obs", "_perm")
                      )
    
    source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
    
    my_x <- merge_dt$tad_rank_obs
    my_y <- merge_dt$tad_rank_perm
    
    outFile <- file.path(outFolder, paste0(dirname(plot_ds), "_", basename(plot_ds), "_TADrank_", rd_type, "_vs_obs_densplot.", plotType))    
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    
    
    densplot(
      x=my_x,
      xlab = "TAD rank: obs.",
      ylab = paste0("TAD rank: ", rd_type),
      main = paste0(plot_ds),
      y=my_y,
      pch=16,
      cex=0.7
    )
    mtext(side=3, text=paste0("with perm. ", rd_type))
    addCorr(x=my_x, y=my_y, legPos = "topleft", bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))  
    
        
    
    
  }
  
  
  
}

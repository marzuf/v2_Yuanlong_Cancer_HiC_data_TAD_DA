
# Rscript check_pseudosignifTADs.R

require(foreach)
require(doMC)
registerDoMC(40)

outFolder <- "CHECK_PSEUDOSIGNIFTADS"
dir.create(outFolder, recursive = TRUE)

nPerm <- 1000

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

buildData <- FALSE

if(buildData) {
  
  check_obs_shuff_dt <- foreach(i_perm = 1:nPerm, .combine='rbind') %dopar% {
    
    dt <- get(load(file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_SHUFFLE", i_perm, "all_data_list.Rdata")))  
    all_ds <- names(dt)
    ds = all_ds[1]
    
    ds_dt <- foreach(ds = all_ds, .combine='rbind') %do% {
      
      sub_dt <- final_dt[final_dt$hicds == dirname(ds) & final_dt$exprds == basename(ds),]
      
      pseudo_signif_tads_dt <- dt[[paste0(ds)]]
      
      stopifnot(setequal(unique(pseudo_signif_tads_dt[["dataset_g2t_dt"]]$region), 
                         unique(pseudo_signif_tads_dt[["dataset_tadpos_dt"]]$region)))
      
      stopifnot(setequal(unique(pseudo_signif_tads_dt[["dataset_geneList"]]),pseudo_signif_tads_dt[["dataset_g2t_dt"]]$entrezID))
      
      pseudo_tads_g2t_dt <- unique(pseudo_signif_tads_dt[["dataset_g2t_dt"]])
      
      pseudo_tads <- unique(pseudo_signif_tads_dt[["dataset_tadpos_dt"]]$region)
      
      pseudo_genes <- unique(pseudo_signif_tads_dt[["dataset_geneList"]])
      
      true_signif_dt <- sub_dt[sub_dt$adjPvalComb <= 0.01,]
      
      g2t_dt <- read.delim(file.path(dirname(ds), "genes2tad", "all_genes_positions.txt"), stringsAsFactors = FALSE,
                           header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))
      
      geneList <- get(load(file.path("PIPELINE/OUTPUT_FOLDER/", ds, "0_prepGeneData", "pipeline_geneList.Rdata")))
      
      stopifnot(geneList %in% g2t_dt$entrezID)
      
      pip_g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
      stopifnot(true_signif_dt$region %in% pip_g2t_dt$region)
      
      true_signif_g2t <- pip_g2t_dt[pip_g2t_dt$region %in% true_signif_dt$region,]
      
      
      stopifnot(pseudo_genes %in% geneList)
      stopifnot(pseudo_genes %in% pip_g2t_dt$entrezID)
      
      stopifnot(length(unique(true_signif_dt$region)) == length(pseudo_tads))
      stopifnot(length(unique(true_signif_g2t$region)) == length(pseudo_tads))
      
      stopifnot(pseudo_tads_g2t_dt$entrezID %in% geneList)
      stopifnot(pseudo_tads_g2t_dt$region %in% pseudo_tads)
      
      pseudo_count <- table(pseudo_tads_g2t_dt$region)
      
      obs_count <- table(true_signif_g2t$region)
      
      stopifnot(length(pseudo_count) == length(obs_count))
      
      rbind(data.frame(
        dataset=ds,
        nGenesInTAD=as.numeric(pseudo_count),
        dataTye="shuffle",
        stringsAsFactors = FALSE
      ),
      data.frame(
        dataset=ds,
        nGenesInTAD=as.numeric(obs_count),
        dataTye="observed",
        stringsAsFactors = FALSE
      ))
    }
    ds_dt
  }
  
  outFile <- file.path(outFolder, "check_obs_shuff_dt.Rdata")
  save(check_obs_shuff_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
} else {
  outFile <- file.path(outFolder, "check_obs_shuff_dt.Rdata")
  check_obs_shuff_dt <- get(load(outFile))
  
}

# load("CHECK_PSEUDOSIGNIFTADS/check_obs_shuff_dt.Rdata")
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

totObs <- sum(check_obs_shuff_dt$nGenesInTAD[check_obs_shuff_dt$dataTye == "observed"])
totShuff <- sum(check_obs_shuff_dt$nGenesInTAD[check_obs_shuff_dt$dataTye == "shuffle"])

png(file.path(outFolder, "nGenes_in_TADs.png"), height=400, width=450)
plot_multiDens(
 list(obs= check_obs_shuff_dt$nGenesInTAD[check_obs_shuff_dt$dataTye == "observed"],
 shuffle= check_obs_shuff_dt$nGenesInTAD[check_obs_shuff_dt$dataTye == "shuffle"]),
 plotTit = "# of genes in signif. TADs"
)
mtext(side=3, text=paste0("tot. # genes obs. = ", totObs, "; tot. # genes shuff. = ", totShuff))
dev.off()

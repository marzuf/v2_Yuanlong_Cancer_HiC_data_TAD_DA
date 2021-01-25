hicds2tissue <- c(
  "Barutcu_MCF-10A_40kb"="Breast_Mammary_Tissue",
  "Barutcu_MCF-7_40kb"="Breast_Mammary_Tissue",
  "ENCSR079VIJ_G401_40kb"="Kidney_Cortex", 
  "ENCSR312KHQ_SK-MEL-5_40kb"="Skin_Not_Sun_Exposed_Suprapubic",
  "ENCSR346DCU_LNCaP_40kb"="Prostate",
  "ENCSR401TBQ_Caki2_40kb"="Kidney_Cortex",  
  "ENCSR444WCZ_A549_40kb"="Lung",
  "ENCSR489OCU_NCI-H460_40kb"="Lung",    
  "ENCSR504OTV_transverse_colon_40kb"="Colon_Transverse", 
  "ENCSR549MGQ_T47D_40kb"="Breast_Mammary_Tissue",        
  "ENCSR862OGI_RPMI-7951_40kb"="Skin_Not_Sun_Exposed_Suprapubic",
  "GSE105194_cerebellum_40kb"="Brain_Cerebellum",
  "GSE105194_spinal_cord_40kb"="Brain_Cerebellum",
  "GSE105318_DLD1_40kb"="Colon_Transverse",  
  "GSE105381_HepG2_40kb"="Liver",   
  "GSE109229_BT474_40kb" ="Breast_Mammary_Tissue",               
  "GSE109229_SKBR3_40kb"="Breast_Mammary_Tissue",   
  "GSE118514_22Rv1_40kb"  ="Prostate",           
  "GSE118514_RWPE1_40kb"="Prostate",
  "GSE118588_Panc_beta_40kb"= "Pancreas",
  "GSE99051_786_O_40kb"="Kidney_Cortex",  
  "HMEC_40kb"="Breast_Mammary_Tissue",      
  "K562_40kb"="Whole_Blood",
  "LG1_40kb"="Lung",       
  "LG2_40kb"="Lung",       
  "LI_40kb"="Liver",           
  "PA2_40kb" = "Pancreas",
  "PA3_40kb" ="Pancreas",    
  "Panc1_rep12_40kb"="Pancreas",
  "Rao_HCT-116_2017_40kb"="Colon_Transverse")
# Rscript cops_signif.R

setDir <- "/media/electron"
setDir <- ""

allgenes_dt <- read.delim(file.path(setDir,
                                    "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                          header=T, stringsAsFactors = FALSE)
allgenes_dt$entrezID <- as.character(allgenes_dt$entrezID)
allgenes_dt$symbol <- as.character(allgenes_dt$symbol)
stopifnot(!duplicated(allgenes_dt$entrezID))
stopifnot(!duplicated(allgenes_dt$symbol))
entrez2symb <- setNames(allgenes_dt$symbol, allgenes_dt$entrezID)

signifThresh <- 0.01

plotType <- "png"
ggHeight <- 5
ggWidth <- 5

cops_dt <- read.delim("all_cops_gtex.tsv", stringsAsFactors = FALSE, header=TRUE)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(ggpubr)
require(foreach)
require(doMC)
registerDoMC(40)

outFolder <- "COPS_SIGNIF"
dir.create(outFolder, recursive = TRUE)

buildTable <- FALSE

signif_dt <- get(load(file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")))
signif_dt$signif_lab <- signif_dt$adjPvalComb <= signifThresh

signif_dt$ds <- file.path(signif_dt$hicds, signif_dt$exprds)
ds = unique(signif_dt$ds)[1]

if(buildTable){
  all_ds_corr_signif_dt <- foreach(ds = unique(signif_dt$ds), .combine='rbind') %dopar% {
    
    sub_signif_dt <- signif_dt[signif_dt$ds == ds,]
    stopifnot(nrow(sub_signif_dt) > 0)
    
    hicds <-  dirname(ds)
    exprds <- basename(ds)
    tissue <- hicds2tissue[hicds]
    stopifnot(!is.na(tissue))
    
    stopifnot(hicds %in% signif_dt$hicds)
    stopifnot(exprds %in% signif_dt$exprds)
    
    ds_signif_dt <- signif_dt[signif_dt$hicds == hicds & signif_dt$exprds == exprds,]
    stopifnot(nrow(ds_signif_dt) > 0)
    region2signif <- setNames(ds_signif_dt$signif_lab, ds_signif_dt$region)
    
    stopifnot(tissue %in% cops_dt$tissue)
    tissue_cops_dt <- cops_dt[cops_dt$tissue == tissue,]
    stopifnot(nrow(tissue_cops_dt) > 0)
    nrow(tissue_cops_dt)
    
    g2t_dt <- read.delim(file.path(hicds,"genes2tad", "all_genes_positions.txt"), header=F,
                         col.names=c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    g2t_dt$symbol <- entrez2symb[g2t_dt$entrezID]
    stopifnot(!is.na(g2t_dt$symbol))
    symbol2region <- setNames(g2t_dt$region, g2t_dt$symbol)
    
    ds_cops_dt <- tissue_cops_dt[tissue_cops_dt$gene1_name %in% g2t_dt$symbol &
                                   tissue_cops_dt$gene2_name %in% g2t_dt$symbol,
    ]
    nrow(ds_cops_dt)
    stopifnot(nrow(ds_cops_dt) > 0)
    
    ds_cops_dt$gene1_region <- symbol2region[ds_cops_dt$gene1_name]
    ds_cops_dt$gene2_region <- symbol2region[ds_cops_dt$gene2_name]
    stopifnot(!is.na(ds_cops_dt))
    
    ds_cops_dt$gene1_region_signifLab <- region2signif[ds_cops_dt$gene1_region]
    ds_cops_dt$gene2_region_signifLab <- region2signif[ds_cops_dt$gene2_region]
    # stopifnot(!is.na(ds_cops_dt)) # TADs that have not been included in the pipeline will be NA
    ds_cops_dt$sameTAD <- ds_cops_dt$gene1_region == ds_cops_dt$gene2_region
    ds_cops_dt$ds <- ds
    ds_cops_dt
  }
  outFile <- file.path(outFolder, "all_ds_corr_signif_dt.Rdata")
  save(all_ds_corr_signif_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "all_ds_corr_signif_dt.Rdata")
  all_ds_corr_signif_dt <- get(load(outFile))
  # load("COPS_SIGNIF/all_ds_corr_signif_dt.Rdata")
}
stopifnot(setequal(all_ds_corr_signif_dt$ds, signif_dt$ds))

nDS <- length(unique(all_ds_corr_signif_dt$ds))

plotTit <- paste("gene pair correlations")

# plot_multiDens(split(all_ds_corr_signif_dt$correlation, all_ds_corr_signif_dt$sameTAD))

outFile <- file.path(outFolder, paste0("corrValues_allDS_sameTAD_vs_diffTAD_boxplot.", plotType))
p1 <- ggboxplot(data = all_ds_corr_signif_dt,
                add = "jitter",
          x = "sameTAD",
          y="correlation") + 
  labs(x="sameTAD", y="correlation")+ 
  ggtitle(plotTit, subtitle=paste0("all pairs - all DS (n=", nDS, ")"))+
 stat_compare_means()
ggsave(p1, filename=outFile, width=ggWidth, height = ggHeight)
cat(paste0("... written: ", outFile, "\n"))


signif_sameTADs_cops_dt <- all_ds_corr_signif_dt[!is.na(all_ds_corr_signif_dt$gene1_region_signifLab) &
                                        !is.na(all_ds_corr_signif_dt$gene2_region_signifLab) &
                                          all_ds_corr_signif_dt$gene2_region == all_ds_corr_signif_dt$gene1_region,
]
stopifnot(signif_sameTADs_cops_dt$gene1_region_signifLab == signif_sameTADs_cops_dt$gene2_region_signifLab)

outFile <- file.path(outFolder, paste0("corrValues_allDS_sameTAD_signif_vs_notsignif_boxplot.", plotType))
p2 <- ggboxplot(data = signif_sameTADs_cops_dt,
                add = "jitter",
          x = "gene1_region_signifLab",
          y="correlation") + 
  labs(x="signif", y="correlation")+ 
  ggtitle(plotTit, subtitle=paste0("pairs from same TADs - all DS (n=", nDS, ")"))+
  stat_compare_means()
ggsave(p2, filename=outFile, width=ggWidth, height = ggHeight)
cat(paste0("... written: ", outFile, "\n"))
  
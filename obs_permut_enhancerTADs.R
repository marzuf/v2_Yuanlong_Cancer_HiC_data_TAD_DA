require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

startTime <- Sys.time()

script_name <- "obs_permut_enhancerTADs.R"

# Rscript obs_permut_enhancerTADs.R

plotType <- "svg"
myHeight <- myWidth <- 7
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "OBS_PERMUT_ENHANCERTADS"
dir.create(outFolder, recursive = TRUE)


final_DT <- get(load(file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")))
final_DT_permut <- get(load(file.path("CREATE_FINAL_TABLE_RANDOM", "all_result_dt.Rdata")))

exprds <- "TCGAluad_norm_luad"
famType <- "hgnc_family_short"

tad_signif_thresh <- 0.05

myhicds <- "ENCSR489OCU_NCI-H460"

all_hicds <- list.files(file.path("PIPELINE/OUTPUT_FOLDER"))

myhicds <- "ENCSR489OCU_NCI-H460"

exprds <- "TCGAluad_norm_luad"

buildData <- TRUE

all_hicds <- all_hicds[grep(myhicds, all_hicds)]

enh2entrez_dt <- get(load("enhanceratlas_data/PREP_TARGETGENES/enhancerEntrezTarget_dt.Rdata"))

if(buildData) {
  
  all_data_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    
    cat(paste0("... start ", hicds, "\n"))
    
    enh2entrez_g2t_dt <- enh2entrez_dt
    colnames(enh2entrez_g2t_dt)[colnames(enh2entrez_g2t_dt) == "target_entrezID"] <- "entrezID"
    enh2entrez_g2t_dt$entrezID <- as.character(enh2entrez_g2t_dt$entrezID)
    
    
    g2t_dt_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    
    sum(enh2entrez_g2t_dt$entrezID %in% g2t_dt$entrezID)/nrow(enh2entrez_g2t_dt)
    
    
    enh2entrez_g2t_dt <- merge(enh2entrez_g2t_dt, g2t_dt[, c("entrezID", "region")], all.x=FALSE, all.y=FALSE, by="entrezID")
    colnames(enh2entrez_g2t_dt)[colnames(enh2entrez_g2t_dt) == "region"] <- "target_region"
    stopifnot(!is.na(enh2entrez_g2t_dt))
    cat(paste0(nrow(enh2entrez_dt) , "->", nrow(enh2entrez_g2t_dt), "\n"))
    
    enh2tad_dt <- get(load(file.path("enhanceratlas_data/PREP_ENHANCER2TAD_EP", paste0(hicds, "_enhancer2tad_dt.Rdata"))))
    colnames(enh2tad_dt)[colnames(enh2tad_dt) == "region"] <- "enhancer_region"
    
    
    enh2tad_dt$enhancer_id <- file.path(enh2tad_dt$enhancer_chromo, enh2tad_dt$enhancer_start,enh2tad_dt$enhancer_end)
    enh2entrez_g2t_dt$enhancer_id <- file.path(enh2entrez_g2t_dt$enhancer_chromo, enh2entrez_g2t_dt$enhancer_start,enh2entrez_g2t_dt$enhancer_end)
    stopifnot(enh2entrez_g2t_dt$enhancer_id %in% enh2tad_dt$enhancer_id)
    
    enh2tad_g2t_dt <- merge(enh2entrez_g2t_dt, unique(na.omit(enh2tad_dt[,c("enhancer_id", "enhancer_region")])), by=c("enhancer_id"))
    
    stopifnot(nrow(enh2tad_g2t_dt) == nrow(enh2entrez_g2t_dt))
    
    stopifnot(!is.na(enh2tad_g2t_dt))
    
    enh2tad_g2t_dt$target_entrez_sameTAD <- enh2tad_g2t_dt$enhancer_region == enh2tad_g2t_dt$target_region
    
    tad_nEPsameTAD_dt <- aggregate(target_entrez_sameTAD ~ target_region, data=enh2tad_g2t_dt, FUN=function(x)sum(x))
    colnames(tad_nEPsameTAD_dt) [colnames(tad_nEPsameTAD_dt) == "target_entrez_sameTAD"] <- "nEPsameTAD"
    tad_nEPdiffTAD_dt <- aggregate(target_entrez_sameTAD ~ target_region, data=enh2tad_g2t_dt, FUN=function(x)sum(!x))
    colnames(tad_nEPdiffTAD_dt) [colnames(tad_nEPdiffTAD_dt) == "target_entrez_sameTAD"] <- "nEPdiffTAD"
    tad_nEP_dt <- aggregate(target_entrez_sameTAD ~ target_region, data=enh2tad_g2t_dt, FUN=length)
    colnames(tad_nEP_dt) [colnames(tad_nEP_dt) == "target_entrez_sameTAD"] <- "nEP"
    
    tad_nEP_all_dt <- merge(tad_nEP_dt, merge(tad_nEPdiffTAD_dt, tad_nEPsameTAD_dt, by="target_region"), by="target_region")
    stopifnot(nrow(tad_nEP_all_dt) == length(unique(enh2tad_g2t_dt$target_region)))
    stopifnot(!is.na(tad_nEP_all_dt))
    
    colnames(tad_nEP_all_dt)[colnames(tad_nEP_all_dt)=="target_region"] <- "region"
    
    if(hicds %in% final_DT$hicds) {
      result_dt <- final_DT[final_DT$hicds == hicds & final_DT$exprds == exprds,]
    } else if(hicds %in% final_DT_permut$hicds) {
      result_dt <- final_DT_permut[final_DT_permut$hicds == hicds & final_DT_permut$exprds == exprds,]
    } else {
      stop("error")
    }
    
    tad_nEP_all_dt_final <- merge(tad_nEP_all_dt, result_dt[,c("hicds", "exprds", "region", "meanCorr", "adjPvalComb")],by="region", all.x=FALSE, all.y=FALSE)
    stopifnot(!is.na(tad_nEP_all_dt_final))
    
    cat(paste0(nrow(result_dt) , " -> ", nrow(tad_nEP_all_dt_final), "\n"))
    tad_nEP_all_dt_final
    
    
  }
  outFile <- file.path(outFolder, "all_data_dt.Rdata")
  save(all_data_dt, file =outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
} else{
  
  outFile <- file.path(outFolder, "all_data_dt.Rdata")
  all_data_dt <- get(load(outFile))
}





### => BY TAD

all_data_dt$adjPvalComb_log10 <- -log10(all_data_dt$adjPvalComb )
all_data_dt$hicds_lab <- gsub("ENCSR489OCU_NCI-H460_(.+)","\\1",  all_data_dt$hicds)

all_data_dt$signif <- ifelse(all_data_dt$adjPvalComb <= tad_signif_thresh, "signif.", "not signif.")


all_data_dt$nEPsameTAD_log10 <- log10(all_data_dt$nEPsameTAD )
all_data_dt$nEPdiffTAD_log10 <- log10(all_data_dt$nEPdiffTAD )

outFile <- file.path(outFolder, paste0("all_", myhicds, "_signifNotSignif_nEPsameTAD_boxplot.", plotType))
p_box <- ggboxplot(all_data_dt, x = "signif", y = "nEPsameTAD_log10") + 
  facet_grid(~hicds_lab, switch="x") 
ggsave(p_box, filename = outFile, height=myHeightGG, width=myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("all_", myhicds, "_signifNotSignif_nEPdiffTAD_boxplot.", plotType))
p_box <- ggboxplot(all_data_dt, x = "signif", y = "nEPdiffTAD_log10") + 
  facet_grid(~hicds_lab, switch="x") 
ggsave(p_box, filename = outFile, height=myHeightGG, width=myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))



for(plot_var in c("nEPdiffTAD", "nEPdiffTAD")) {
  
  for(hicds in all_hicds){
    plot_dt <- all_data_dt[all_data_dt$hicds == hicds,]
    
    
    
    my_x <- plot_dt[,c(plot_var)]
    my_y <- plot_dt[,c("adjPvalComb_log10")]
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", plot_var, "_vs_adjPvalComb_log10_densplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    
    
    densplot(
      x = my_x,
      y=my_y,
      xlab = paste0(plot_var),
      ylab = paste0("TAD adjPvalComb [-log10]"),  
      cex=0.7,
      main = paste0(hicds, "  - ",exprds)
    )  
    addCorr(x=my_x, y=my_y, bty="n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
  }

outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_allTADs_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(all_data_dt[,paste0(plot_var)], all_data_dt$hicds_lab),
               plotTit=plot_var, legPos = "topleft")
mtext(side=3, text = paste0(myhicds, " -", exprds))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_allTADs_boxplot.", plotType))

p_box <- ggboxplot(data=all_data_dt, x="hicds_lab", y=paste0(plot_var))
ggsave(p_box, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

}




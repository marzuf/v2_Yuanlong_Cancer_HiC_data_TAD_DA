
require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)

# Rscript nsignif_and_auccoexprdist_obs_permut.R

plotType <- "png"
myHeight <- myWidth <- 400

plotType <- "svg"
myHeight <- myWidth <- 7

plotTypeGG <- "svg"
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "NSIGNIF_AND_AUCCOEXPRDIST_OBS_PERMUT"
dir.create(outFolder, recursive = TRUE)

tadSignifThresh <- 0.01

result_dt1 <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
result_dt2 <- get(load("CREATE_FINAL_TABLE_RANDOM//all_result_dt.Rdata"))
result_dt <- rbind(result_dt1, result_dt2)

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
all_hicds <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
hicds = all_hicds[1]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))


famType <- "hgnc"
famTypeSub <- "hgnc_family_short"

all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    nSignif <- sum(result_dt$adjPvalComb[result_dt$hicds == hicds & result_dt$exprds == exprds] <= tadSignifThresh)
    
    auc_sameTAD <- get(load(file.path("AUC_COEXPRDIST_WITHFAM_SORTNODUP", hicds, paste0(exprds, "_", famType), famTypeSub, "auc_sameTAD_distVect.Rdata")))
    auc_diffTAD <- get(load(file.path("AUC_COEXPRDIST_WITHFAM_SORTNODUP", hicds, paste0(exprds, "_", famType), famTypeSub, "auc_diffTAD_distVect.Rdata")))
    all_auc_values <- get(load(file.path("AUC_COEXPRDIST_WITHFAM_SORTNODUP", hicds, paste0(exprds, "_", famType), famTypeSub, "auc_values.Rdata")))
    
    stopifnot(all_auc_values[["auc_sameTAD_distVect"]]/all_auc_values[["auc_diffTAD_distVect"]] == auc_sameTAD/auc_diffTAD)
    stopifnot(all_auc_values[["auc_sameTAD_distVect"]]/all_auc_values[["auc_diffTAD_distVect"]] == all_auc_values[["auc_ratio_same_over_diff_distVect"]])
    
    data.frame(
      hicds=hicds,
      exprds=exprds,
      nSignifTADs=nSignif,
      auc_ratio_sameTAD_diffTAD = auc_sameTAD/auc_diffTAD,
      stringsAsFactors = FALSE
    )
      }
  exprds_dt
}

rd_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMSHIFT", "PERMUTG2T", "RANDOMMIDPOSDISC" )

all_dt$hicds_sub <- gsub("_RANDOM.+_40kb", "_40kb", all_dt$hicds)
all_dt$hicds_sub <- gsub("_PERMUT.+_40kb", "_40kb", all_dt$hicds)

all_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                          gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                               gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC",
                                    gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                         gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_dt$hicds)))))
all_dt$hicds_lab[! all_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"
all_dt$dataset_sub <- file.path(all_dt$hicds_sub, all_dt$exprds)

outFile <- file.path(outFolder, "all_dt.Rdata")
save(all_dt, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

require(ggpubr)
require(ggsci)

tmp_dt <- all_dt[all_dt$hicds_lab == "OBSERVED",]
tmp_dt <- tmp_dt[order(tmp_dt$nSignifTADs, decreasing = TRUE),]
all_dt$dataset_sub <- factor(all_dt$dataset_sub, levels = tmp_dt$dataset_sub)

nDS <- length(unique(all_dt$dataset_sub))

bar_plot <- ggbarplot(data = all_dt, 
                      y = "nSignifTADs", 
                      x = "hicds_lab" , 
                      fill = "hicds_lab", 
                      title ="# signif. TADs",
                      position=position_dodge(),
                      # xlab = "", 
                      xlab = "Datasets ranked by decreasing obs. # of signif. TADs",
                      ylab="# signif. TADs") + 
  labs(subtitle = paste0("all DS - n = ", nDS, "; TAD signif. p-val <= ", tadSignifThresh), fill ="")+
  scale_fill_nejm()+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  facet_grid(~dataset_sub, switch = "x") +
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"))

outFile <- file.path(outFolder, paste0("cmp_signif_genes_obs_permut_barplot.", plotTypeGG))
ggsave(bar_plot, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


all_types <- unique(all_dt$hicds_lab) 
rd_type = all_types[1]
for(rd_type in all_types) {
  
  plot_dt <- all_dt[all_dt$hicds_lab == rd_type,]
  
  
  
  par(bty ="L")
  plot(
    x = plot_dt$auc_ratio_sameTAD_diffTAD,
    y = plot_dt$nSignifTADs,
    xlab  = "AUC ratio sameTAD/diffTAD",
    ylab ="# signif TADs",
    pch = 16,
    cex = 0.7,
    cex.main = plotCex,
    cex.axis = plotCex,
    cex.lab =plotCex
  )
  
}

 


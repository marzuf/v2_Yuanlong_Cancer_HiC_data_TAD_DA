########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript conservTADs_and_purity.R\n"))

script_name <- "conservTADs_and_purity.R"

# do the same as for GIMAPs but for all conserved TADs !!!


####### !!!!!!!!!!!  FOR THE MOMENT I DO NOT RESTRICT 
# the datasets in which I look at expression (I do not ensure that is a DS where the conserved region is signif. DA)
# rationale: I want to look if this set of genes is related to purity, irrespective of DA


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

# Rscript conservTADs_and_purity.R
# Rscript conservTADs_and_purity.R EPIC

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- TRUE

myHeight <- 400 
myWidth <- 400
plotType <- "png"
plotCex <- 1.4
mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")


plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 11
plotCex <- 1.2



args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  purity_ds <- args[1]  
} else{
  purity_ds <- ""
}
outFolder <- file.path("CONSERVTADS_AND_PURITY", purity_ds)
dir.create(outFolder, recursive = TRUE)

# SET HERE ALL PURITY FILE
all_ds_corrPurity_dt <- get(load(file.path("ALLTADS_AND_PURITY", purity_ds, "all_ds_corrPurity_dt.Rdata")))
all_ds_corrPurity_dt$regID <- file.path(all_ds_corrPurity_dt$dataset, all_ds_corrPurity_dt$region)


# HARD CODED SETTING INPUT FILE
# take only the genes that are at the intersect
all_intersect_genes <- get(load(file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2",
                              "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))

all_conserv <- get(load(file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2",
                                          "conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))

nConserv_dt <- data.frame(conserved_region = names(all_conserv),
                          nConserv = as.numeric(lengths(all_conserv)), stringsAsFactors = FALSE)


all_intersect_genes_dt <- foreach(i = 1:nrow(all_intersect_genes), .combine='rbind') %dopar% {
  entrezIDs <- unlist(strsplit(all_intersect_genes$intersect_genes_entrez[i], split = ","))
  regIDs <- unlist(strsplit(all_intersect_genes$corresp_tads[i], split = ","))
  out_dt <- data.frame(
    conserved_region = all_intersect_genes$conserved_region[i],
    entrezID = rep(entrezIDs, length(regIDs)),
    regID = rep(regIDs, each=length(entrezIDs)),
    stringsAsFactors = FALSE
  )
  stopifnot(!duplicated(file.path(out_dt$regID, out_dt$entrezID)))
  out_dt
}

all_ds_corrPurity_dt$conserved <- all_ds_corrPurity_dt$regID %in% all_intersect_genes_dt$regID


boxplot(all_ds_corrPurity_dt$purityCorr~all_ds_corrPurity_dt$conserved)


all_dt <- merge(all_intersect_genes_dt, all_ds_corrPurity_dt[,c("regID", "entrezID", "purityCorr")],
                by=c("regID", "entrezID"), 
                # all.x=FALSE, all.y=FALSE)
                all.x=TRUE, all.y=FALSE)

# not sure
# stopifnot(!is.na(all_dt))
# > all_dt[28,]
# regID entrezID               region purityCorr
# 28 K562_40kb/TCGAlaml_wt_mutFLT3/chr1_TAD75      712 conserved_region_142         NA
# will need to check with the plot color-code gene conservation

all_dt <- merge(all_dt, nConserv_dt, by="conserved_region", all=TRUE)
# stopifnot(!is.na(all_dt))
# which(apply(all_dt, 1, function(x) any(is.na(x))))


outFile <- file.path(outFolder, "all_dt.Rdata")
save(all_dt, file = outFile, version=2)

plotTit <- "correlation expression of genes at the intersect"

# boxplot of something like
plot(all_dt$purityCorr~all_dt$nConserv)

all_dt <- all_dt[order(all_dt$nConserv),]
all_dt$conserved_region <- factor(all_dt$conserved_region, levels=unique(all_dt$conserved_region))

subTit <- paste0("")


barbox_p <- ggplot(data=all_dt, 
                   aes(x=conserved_region, y =  purityCorr))+ 
  ggtitle(plotTit, subtitle = subTit)+
  geom_boxplot()
  
  
  geom_bar(stat="identity", alpha=0.2) + 
  labs(fill="", color="", x="Conserved in \"x\" datasets", y="# of regions")+
  scale_color_manual(values=mycols, breaks = names(mycols), labels=names(mycols), drop=F)+
  # scale_y_continuous(expand=c(0,5))+
  scale_fill_manual(values=mycols, drop=F)+
  geom_boxplot(data=nConsByPermut_dt, aes(x=conserved, y = nRegions), 
               color = mycols["shuffle"],
               inherit.aes = FALSE, notch = TRUE)+
  theme(
    text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.line = element_line(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold"),
    legend.text = element_text(size=12)
  ) 

outFile <- file.path(outFolder, paste0("obsCons_permCons_nCons_barplot_boxplot.", plotType))
ggsave(barbox_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))





stop("--ok<")



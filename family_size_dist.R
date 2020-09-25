

# Rscript family_size_dist.R

# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity

set.seed(20200903) # NB forgot the first time I launched it

require(doMC)
require(foreach)
registerDoMC(40)
require(reshape2)
require(igraph)
require(ggpubr)

runFolder <- "."
pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")
familyVar <- "hgnc_family_short"

corMethod <- "pearson"
familyData <- "hgnc_family_short"

nRandom <- 100

setDir <- "/media/electron"
setDir <- ""


aggFunc <- "mean"

plotType <- "svg"
myHeight <- 5
myWidth <- 7

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- file.path("FAMILY_SIZE_DIST")
# outFolder <- file.path("SAMEFAMINSAMETAD")
dir.create(outFolder, recursive = TRUE)

hgnc_geneFamilyFile <- file.path(setDir, "/mnt/ed4/marie/family_data_2/hgnc_entrez_family.txt")
hgnc_geneFamilyDT <- read.delim(hgnc_geneFamilyFile, col.names=c("entrezID", "family"), header = F, stringsAsFactors = F)
hgnc_geneFamilyDT$entrezID <- as.character(hgnc_geneFamilyDT$entrezID)
hgnc_geneFamilyDT$family_short <- unlist(sapply(hgnc_geneFamilyDT$family, function(x) strsplit(x, "\\|")[[1]][1] ))
# any(duplicated(hgnc_geneFamilyDT$entrezID))
# FALSE
familyDist_dt <- data.frame(
  family=as.character(names(table(hgnc_geneFamilyDT$family_short))),
  nGenes=as.numeric(table(hgnc_geneFamilyDT$family_short)),
  stringsAsFactors = FALSE
)
familyDist_dt$nGenes_log10 <- log10(familyDist_dt$nGenes)
familyDist_dt <- familyDist_dt[order(familyDist_dt$nGenes, decreasing = TRUE),]

head(familyDist_dt)

plotTit <- "# genes by family"
mySub <- paste0("(# families = ", nrow(familyDist_dt), ")")

p2 <- ggdensity(familyDist_dt,
                fill = "lightgray",
                x = "nGenes_log10",
                y = "..density..",
                # combine = TRUE,                  # Combine the 3 plots
                xlab = "# genes by family [log10]", 
                # add = "median",                  # Add median line. 
                rug = TRUE,                      # Add marginal rug
                # color = "signif",
                # fill = "signif",
                palette = "jco"
) +
  # scale_color_manual(values=my_cols)+
  # scale_fill_manual(values=my_cols)  +
  ggtitle(plotTit, subtitle = mySub)+
  # labs(color=paste0(legTitle),fill=paste0(legTitle)
  labs( y="Density") + 
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold")
  ) 

outFile <- file.path(outFolder, paste0("nbrGenes_by_family_log10.", plotType))
ggsave(p2, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("familyDist_dt.Rdata"))
save(familyDist_dt, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


curr_var <- "nbrEdges"
inFile <- file.path("SAMEFAMINSAMETAD_V2_ALLDS", "all_ds_results.Rdata")
all_ds_results <- get(load(inFile))
stopifnot(length(all_ds_results) == 58)
all_obs_nbr <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[[paste0(curr_var)]])))

# aggregate by family
nbr_obs_dt <- data.frame(family_lab=names(eval(parse(text=paste0("all_obs_nbr")))), 
                         nbr=as.numeric(eval(parse(text=paste0("all_obs_nbr")))), stringsAsFactors = FALSE)
nbr_obs_dt$family <- gsub(".+/.+?\\.(.+)", "\\1", nbr_obs_dt$family_lab)
stopifnot(nbr_obs_dt$family %in% familyDist_dt$family)
stopifnot(nbr_obs_dt$family != "")
stopifnot(!is.na(nbr_obs_dt$family))
agg_obs_dt <- aggregate(nbr~family, data=nbr_obs_dt, FUN=aggFunc)
stopifnot(agg_obs_dt$family %in% familyDist_dt$family)
colnames(agg_obs_dt)[colnames(agg_obs_dt) == "nbr"] <- paste0(curr_var)

logOffset <- 0.1

plotCex <- 1.2

agg_obs_dt[,paste0(curr_var, "_log10")] <- log10(logOffset + agg_obs_dt[,paste0(curr_var)])

plot_dt <- merge(agg_obs_dt, familyDist_dt, all.x=T, all.y=F, by="family")
stopifnot(!is.na(plot_dt))

# plot_dt <- merge(agg_obs_dt, agg_rd_dt, by="family", all=T, suffixes=c("_obs", "_rd"))
# stopifnot(!is.na(plot_dt))

plotTit <- paste0(curr_var, " and family size")
subTit <- paste0(curr_var, " aggreg. func = ", aggFunc, "; # families = ", length(unique(plot_dt$family)))

my_x <- plot_dt[,paste0( "nGenes_log10")]
my_y <- plot_dt[,paste0(curr_var, "_log10")]

outFile <- file.path(outFolder,  paste0("allDS_", curr_var, "_byFam_vs_nbrGenes_scatterplot_", gsub("\\.", "", logOffset), ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight))
plot(x=my_x,y=my_y, main=plotTit, 
     pch=16, cex=0.7,
     xlab=paste0("# genes [log10]"), 
     cex.axis=plotCex,
     cex.lab=plotCex,
     cex.main=plotCex,
     ylab=paste0("# ", curr_var, " [log10(+", logOffset, ")]"))
# abline(lm(my_y~my_x))
addCorr(x=my_x, y=my_y, legPos="topleft", bty="n")
mtext(side=3, text = subTit, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))





# nbr_rd_dt <- data.frame(family_lab=names(eval(parse(text=paste0("all_random_nbr", curr_var)))), 
#                         nbr=as.numeric(eval(parse(text=paste0("all_random_nbr", curr_var)))), 
#                         stringsAsFactors = FALSE)
# nbr_rd_dt$family <- gsub(".+/.+\\.(.+)\\.result.+", "\\1", nbr_rd_dt$family_lab)
# stopifnot(nbr_rd_dt$family != "")
# stopifnot(!is.na(nbr_rd_dt$family))
# 
# stopifnot(setequal(nbr_obs_dt$family, nbr_rd_dt$family))
# agg_rd_dt <- aggregate(nbr~family, data=nbr_rd_dt, FUN=aggFunc)








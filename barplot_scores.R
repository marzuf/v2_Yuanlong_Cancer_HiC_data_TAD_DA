


options(scipen=100)

# Rscript barplot_scores.R

script_name <- "barplot_scores.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(ggplot2)
require(reshape2)

registerDoMC(4)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myWidthGG <- 12
myHeightGG <- 8

fcc_col <- "dodgerblue3"
coexpr_col <- "goldenrod"

script17_name <- "170revision2EZH2_score_auc_pval_permGenes"

famType1 <- "hgnc"
famType2 <- "hgnc_family_short"
auc_coexprdist_fold <- file.path("AUC_COEXPRDIST_WITHFAM_SORTNODUP") 
stopifnot(dir.exists(auc_coexprdist_fold))


setDir <- "/media/electron"
setDir <- ""

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

outFolder <- "BARPLOT_SCORES"
dir.create(outFolder, recursive = TRUE)

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

cat(paste0("n allDS = ", length(all_datasets), "\n"))

minTADsize <- 3

hicds = all_hicds[1]

all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    fcc_auc_file <- file.path(pipFolder, hicds, exprds, script17_name, "auc_ratios.Rdata")
    stopifnot(file.exists(fcc_auc_file))
    fcc_auc <- as.numeric(get(load(fcc_auc_file))["prodSignedRatio_auc_permGenes"])
    stopifnot(!is.na(fcc_auc))
    coexpr_auc_file <- file.path(auc_coexprdist_fold, hicds, paste0(exprds, "_", famType1), famType2, "auc_values.Rdata")
    stopifnot(file.exists(coexpr_auc_file))
    coexpr_auc <- get(load(coexpr_auc_file))[["auc_ratio_same_over_diff_distVect"]]
    stopifnot(!is.na(coexpr_auc))
    data.frame(hicds=hicds, exprds=exprds,fcc_auc=fcc_auc, coexpr_auc=coexpr_auc, stringsAsFactors = FALSE)
  } # end-foreach iterating over exprds
  exprds_dt
} # end-foreach iterating over hicds

save(all_dt, file=file.path(outFolder, "all_dt.Rdata"),version=2)

all_dt$dataset <- paste0(all_dt$hicds, "\n", all_dt$exprds)
all_dt <- all_dt[order(all_dt$fcc_auc, decreasing=TRUE),]
ds_levels_fcc <- as.character(all_dt$dataset)

all_dt <- all_dt[order(all_dt$coexpr_auc, decreasing=TRUE),]
ds_levels_coexpr <- as.character(all_dt$dataset)

all_dt <- melt(all_dt, id=c("hicds", "exprds", "dataset"))

### => FCC sorted

all_dt$dataset <- factor(all_dt$dataset, levels=ds_levels_fcc)
all_dt <- all_dt[order(as.numeric(all_dt$dataset)),]
all_dt$exprds_type_col <- all_cols[all_cmps[all_dt$exprds]]
mycols <- all_dt$exprds_type_col[all_dt$variable == "fcc_auc"]

p_var <-  ggplot(all_dt, aes(x = dataset, y = value, fill = variable)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE) +
  ggtitle("coexpr. and FCC AUC", subtitle = "(FCC sorted)")+
  scale_x_discrete(name="")+
  labs(fill="")+
  scale_fill_manual(values=c(fcc_auc=fcc_col, coexpr_auc=coexpr_col), labels=c("FCC", "coexpr."))+
  scale_y_continuous(name=paste0("AUC"),
                     breaks = scales::pretty_breaks(n = 10))+
  geom_hline(yintercept=1, linetype=2)+
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    strip.text.x = element_text(size = 10),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .3, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.text.x = element_text(color=mycols, hjust=1,vjust = 0.5, size=7, angle=90),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )

outFile <- file.path(outFolder, paste0("all_ds_auc_scores_barplot_fccSorted.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


### => coexpr sorted

all_dt$dataset <- factor(all_dt$dataset, levels=ds_levels_coexpr)
all_dt <- all_dt[order(as.numeric(all_dt$dataset)),]
all_dt$exprds_type_col <- all_cols[all_cmps[all_dt$exprds]]
mycols <- all_dt$exprds_type_col[all_dt$variable == "coexpr_auc"]

p_var <-  ggplot(all_dt, aes(x = dataset, y = value, fill = variable)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE) +
  ggtitle("coexpr. and FCC AUC", subtitle = "(coexpr. sorted)")+
  scale_x_discrete(name="")+
  labs(fill="")+
  scale_y_continuous(name=paste0("AUC"),
                     breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values=c(fcc_auc=fcc_col, coexpr_auc=coexpr_col), labels=c("FCC", "coexpr."))+
  geom_hline(yintercept=1, linetype=2)+
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    strip.text.x = element_text(size = 10),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .3, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.text.x = element_text(color=mycols, hjust=1,vjust = 0.5, size=7, angle=90),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )

outFile <- file.path(outFolder, paste0("all_ds_auc_scores_barplot_coexprSorted.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))





#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))






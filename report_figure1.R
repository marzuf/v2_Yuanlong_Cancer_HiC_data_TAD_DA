# FIGURE 1: barplot # TADs, # signif. TADs, # signif. TADs and FDR

options(scipen=100)

setDir = ""

# Rscript report_figure1.R   # 

script_name <- "report_figure1.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(ggplot2)
require(ggpubr)
require(ggforce)
require(ggsci)
require(reshape2)
require(foreach)
require(doMC)
require(plotrix)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

all_cols[all_cols == "red"] <- "brown3"
all_cols[all_cols == "blue"] <- "darkblue"
all_cols[all_cols == "green"] <- "forestgreen"

# col1 <- get_palette("Dark2", 4)[1]
# col2 <- get_palette("Dark2", 4)[2]
# col3 <- get_palette("Dark2", 4)[3]
# col4 <- get_palette("Dark2", 4)[4]

col1 <- pal_d3()(4)[1]
col2 <- pal_d3()(4)[2]
col3 <- pal_d3()(4)[3]
col4 <- pal_d3()(4)[4]


plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- ifelse(plotType=="png", 600, 10)
plotCex <- 1.2

myWidthGG <- 12
myHeightGG <- 8

outFolder <- "REPORT_FIGURE1"
dir.create(outFolder, recursive=TRUE)

signifTAD_thresh <- 0.05
signifTAD_FDR <- 0.2

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

final_table_DT$pval_signif <- final_table_DT$adjPvalComb <= signifTAD_thresh
final_table_DT$fdr_signif <- final_table_DT[, paste0("signifFDR_", signifTAD_FDR)]
final_table_DT$pval_fdr_signif <- final_table_DT$fdr_signif & final_table_DT$pval_signif

count_tot_dt <- aggregate(pval_signif ~ hicds + exprds, data=final_table_DT, FUN=length)
colnames(count_tot_dt)[colnames(count_tot_dt) == "pval_signif"] <- "tot"
count_tot_dt_m <- melt(count_tot_dt, id=c("hicds", "exprds"))
count_pval_dt <- aggregate(pval_signif ~ hicds + exprds, data=final_table_DT, FUN=sum)
count_pval_dt_m <- melt(count_pval_dt, id=c("hicds", "exprds"))
count_fdr_dt <- aggregate(fdr_signif ~ hicds + exprds, data=final_table_DT, FUN=sum)
count_fdr_dt_m <- melt(count_fdr_dt, id=c("hicds", "exprds"))
count_pval_fdr_dt <- aggregate(pval_fdr_signif ~ hicds + exprds, data=final_table_DT, FUN=sum)
count_pval_fdr_dt_m <- melt(count_pval_fdr_dt, id=c("hicds", "exprds"))

count_tot_dt <- count_tot_dt[order(count_tot_dt$tot, decreasing = TRUE),]
count_tot_dt$dataset <- paste0(count_tot_dt$hicds, "\n", count_tot_dt$exprds)
ds_levels <- as.character(count_tot_dt$dataset)
stopifnot(!duplicated(ds_levels))

stopifnot(count_tot_dt$exprds %in% names(all_cmps))

ds_cols <- all_cols[all_cmps[paste0(count_tot_dt$exprds)]]
stopifnot(!is.na(ds_cols))


plot_dt <- rbind(count_tot_dt_m, count_pval_dt_m, count_fdr_dt_m, count_pval_fdr_dt_m)
plot_dt$dataset <- paste0(plot_dt$hicds, "\n", plot_dt$exprds)
plot_dt$dataset <- factor(plot_dt$dataset, levels = ds_levels)
stopifnot(!is.na(plot_dt$dataset))

plot_dt$variable <- as.character(plot_dt$variable)

colsOrder_0 <- c("tot", "pval_signif", "fdr_signif", "pval_fdr_signif")
pval_signif_name <- paste0("pval <= ", signifTAD_thresh)
fdr_signif_name <- paste0("FDR <= ", signifTAD_FDR)
colsOrder <- c("tot", 
               pval_signif_name, 
               fdr_signif_name,
               paste0(pval_signif_name, "&\n",fdr_signif_name))
plot_dt$variable[plot_dt$variable == "pval_signif"] <- pval_signif_name
plot_dt$variable[plot_dt$variable == "fdr_signif"] <- fdr_signif_name
plot_dt$variable[plot_dt$variable == "pval_fdr_signif"] <-   paste0(pval_signif_name, "&\n",fdr_signif_name)

plot_dt$variable <- factor(plot_dt$variable, levels=colsOrder)

plot_dt <- plot_dt[order(as.numeric(plot_dt$variable)),]

plot_dt$dataset_circle <- "\u25CF"


nDS <- length(ds_levels)

p_var <-  ggplot(plot_dt, aes(x = dataset, y = value, fill = variable)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE) +
  ggtitle("Total # TADs and # signif. TADs", subtitle = paste0("# datasets = ", nDS, " (p-val. thresh. = ", signifTAD_thresh, ")"))+
  # scale_x_discrete(name="")+
  labs(fill="")+
  
  scale_x_discrete(name="", labels=rep("\u25CF", length(unique(plot_dt$dataset))))+
  
  scale_fill_manual(values=c(col1,col2, col3, col4))+
  scale_y_continuous(name=paste0("# signif. TADs"),
                     breaks = scales::pretty_breaks(n = 10))+
  theme( # Increase size of axis lines
    # strip.text = element_text(size = 12),
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
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=10),
    axis.text.x = element_text(color=ds_cols, hjust=1,vjust = 0.5, size=7, angle=90),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )

outFile <- file.path(outFolder, paste0("nbrTADs_tot_and_signif_ggbar.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

stop("-ok\n")

plot_dt2 <- plot_dt[plot_dt$variable != "tot",]
# ds_levels2 <- plot_dt2$dataset[order(plot_dt2$value[plot_dt2$variable=="pval_signif"])]
ds_levels2 <- plot_dt2$dataset[order(plot_dt2$value[plot_dt2$variable== pval_signif_name], decreasing=TRUE)]
plot_dt2$dataset <- factor(as.character(plot_dt2$dataset), levels=ds_levels2)
plot_dt2 <- plot_dt2[order(as.numeric(plot_dt2$dataset)),]

p_var2 <-  ggplot(plot_dt2, aes(x = dataset, y = value, fill = variable)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE) +
  ggtitle("Total # of TADs and # signif. TADs", subtitle = paste0("# datasets = ", nDS, " (p-val. thresh. = ", signifTAD_thresh, ")"))+
  scale_x_discrete(name="")+
  labs(fill="")+
  scale_fill_manual(values=c(col2, col3,col4))+
  scale_y_continuous(name=paste0("# signif. TADs"),
                     breaks = scales::pretty_breaks(n = 10))+
  theme( # Increase size of axis lines
    # strip.text = element_text(size = 12),
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
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=10),
    axis.text.x = element_text(color=ds_cols, hjust=1,vjust = 0.5, size=7, angle=90),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )
outFile <- file.path(outFolder, paste0("nbrTADs_only_signif_ggbar.", plotType))
ggsave(plot = p_var2, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


#################### barplot with breaks
from_gap <- max(plot_dt$value[plot_dt$variable != "tot"]) + 10
to_gap <- min(plot_dt$value[plot_dt$variable == "tot"]) - 10
stopifnot(to_gap > from_gap)


barplot_dt <- merge(count_pval_dt, merge(count_fdr_dt, 
                                         merge(count_tot_dt, count_pval_fdr_dt, by=c("hicds", "exprds")), 
                                         by=c("hicds", "exprds")),
                                         by=c("hicds", "exprds"))
barplot_dt[,"hicds"] <- NULL
barplot_dt[,"exprds"] <- NULL
barplot_dt <- barplot_dt[order(barplot_dt$tot, decreasing = TRUE),]
rownames(barplot_dt) <- barplot_dt$dataset
barplot_dt$dataset <- NULL
# barplot(t(barplot_dt), beside = TRUE)
barplot_dt <- barplot_dt[,paste0(colsOrder_0)]

bar_plot_dt <- t(barplot_dt)

nGaps <- ncol(bar_plot_dt) * (nrow(bar_plot_dt)+1) - 1

# Hack for grouping (leaves the extra space at the end)
bar_plot_dt_withGroups <- as.vector(rbind(bar_plot_dt, rep(NA, ncol(bar_plot_dt))))[1:nGaps]

barcols <- all_cols[all_cmps[paste0(gsub(".+\n(.+)", "\\1", colnames(bar_plot_dt)))]]
stopifnot(!is.na(barcols))

outFile <- file.path(outFolder, paste0("nbrTADs_tot_and_signif_breakbar.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="n")
a = gap.barplot(ceiling(as.matrix(bar_plot_dt_withGroups)), 
                gap=c(from_gap,to_gap),
                xlab="",
                ylab="# of TADs",
                main = "Total # of TADs and # signif. TADs",
                col=rep(c(col1,col2,col3,col4,1), ncol(bar_plot_dt)),  # add 1 for the gap
                ytics = pretty(as.numeric(bar_plot_dt), n=20),
                xaxt='n') # disable the default x-axis

x_pos <- seq(2.5, by=5, length.out = ncol(bar_plot_dt))
# add axis labels at mean position
# axis(1, at=x_pos, colnames(bar_plot_dt), las=2, cex.axis=0.5, line=-1, lwd=0, lwd.ticks=1)
axis(1, at=x_pos, las=2, cex.axis=0.5, line=-1, lwd=0, lwd.ticks=1, labels=F)
mtext(side=1, at=x_pos, colnames(bar_plot_dt), las=2, col = barcols, cex = 0.4, line=-0.5)

legend("topright", legend = colsOrder,
       bty="n",  
       fill=c(col1,col2,col3,col4)) 
axis.break(2, from_gap, breakcol="snow", style="gap")
axis.break(2, from_gap*(1+0.02), breakcol="black", style="slash")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


bar_plot_dt_noTot <- bar_plot_dt[rownames(bar_plot_dt) != "tot", ,drop=F]
bar_plot_dt_noTot <- bar_plot_dt_noTot[,order(bar_plot_dt_noTot["pval_signif",], decreasing = TRUE)]
# bar_plot_dt_noTot <- bar_plot_dt_noTot[,order(bar_plot_dt_noTot[pval_signif_name,], decreasing = TRUE)]
nGaps_noTot <- ncol(bar_plot_dt_noTot) * (nrow(bar_plot_dt_noTot)+1) - 1
bar_plot_dt_withGroups_noTot <- as.vector(rbind(bar_plot_dt_noTot, rep(NA, ncol(bar_plot_dt_noTot))))[1:nGaps_noTot]

outFile <- file.path(outFolder, paste0("nbrTADs_signifOnly_breakbar.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="n")
# barplot(ceiling(as.matrix(bar_plot_dt_withGroups_noTot)), 
#                 beside = TRUE,
#                 xlab="",
#                 ylab="# of TADs",
#                 main = "Total # of TADs and # signif. TADs",
#                 col=rep(c(col2,col3,col4,1), ncol(bar_plot_dt_noTot)),  # add 1 for the gap
#                 # ytics = pretty(as.numeric(bar_plot_dt), n=20),
#                 xaxt='n') # disable the default x-axis
from_gap <- max(bar_plot_dt_noTot[rownames(bar_plot_dt_noTot) != "pval_signif", ]) + 2
to_gap <- min(bar_plot_dt_noTot[rownames(bar_plot_dt_noTot) == "pval_signif", ]) - 2
stopifnot(to_gap > from_gap)

a = gap.barplot(ceiling(as.matrix(bar_plot_dt_withGroups_noTot)), 
                gap=c(from_gap,to_gap),
                xlab="",
                ylab="# of TADs",
                main = "# signif. TADs",
                col=rep(c(col2,col3,col4,1), ncol(bar_plot_dt_noTot)),  # add 1 for the gap
                ytics = pretty(as.numeric(bar_plot_dt_noTot), n=20),
                xaxt='n') # disable the default x-axis

x_pos <- seq(2.5, by=5, length.out = ncol(bar_plot_dt))
# add axis labels at mean position
# axis(1, at=x_pos, colnames(bar_plot_dt), las=2, cex.axis=0.5, line=-1, lwd=0, lwd.ticks=1)
axis(1, at=x_pos, las=2, cex.axis=0.5, line=-1, lwd=0, lwd.ticks=1, labels=F)
mtext(side=1, at=x_pos, colnames(bar_plot_dt), las=2, col = barcols, cex = 0.4, line=-0.5)

legend("topright", legend = colsOrder[colsOrder!="tot"],
       bty="n",  
       fill=c(col2,col3,col4)) 
axis.break(2, from_gap, breakcol="snow", style="gap")
axis.break(2, from_gap*(1+0.02), breakcol="black", style="slash")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


p_var <-  ggplot(plot_dt, aes(x = dataset, y = value, fill = variable)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE) +
  ggtitle("Total # TADs and # signif. TADs", subtitle = paste0("# datasets = ", nDS, " (p-val. thresh. = ", signifTAD_thresh, ")"))+
  scale_x_discrete(name="")+
  labs(fill="")+
  scale_fill_manual(values=c(col1,col2, col3, col4))+
  scale_y_continuous(name=paste0("# signif. TADs"),
                     breaks = scales::pretty_breaks(n = 10))+
  theme( # Increase size of axis lines
    # strip.text = element_text(size = 12),
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
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=10),
    axis.text.x = element_text(color=ds_cols, hjust=1,vjust = 0.5, size=7, angle=90),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )

p_var <- p_var + facet_zoom(y=variable!="tot", zoom.data = ifelse(variable=="tot", FALSE, NA))

outFile <- file.path(outFolder, paste0("nbrTADs_tot_and_zoomSignif_ggbar.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))



##########################################################################################
##########################################################################################
##########################################################################################

cat("*** DONE - ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))







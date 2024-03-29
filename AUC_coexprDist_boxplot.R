startTime <- Sys.time()
cat(paste0("> Rscript AUC_coexprDist_boxplot.R\n"))

# Rscript AUC_coexprDist_boxplot.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich
# Rscript AUC_coexprDist_boxplot.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad

options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

registerDoMC(40)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

hicds = "ENCSR079VIJ_G401_40kb"
exprds= "TCGAkich_norm_kich"


plotType <- "svg"
myHeightGG <- myWidthGG <- 7

plotType <- "png"
myHeightGG <- myWidthGG <- 7

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
hicds <- args[1]
exprds <- args[2]

outFolder <- file.path("AUC_COEXPRDIST_BOXPLOT", hicds, exprds)
dir.create(outFolder, recursive = TRUE)

familyData <- "hgnc"
familySubData <- "hgnc_family_short"

inFolder <- file.path("AUC_COEXPRDIST_WITHFAM_SORTNODUP")

maxDist <- 500*1000

### HARD CODEDrequire(ggplot2)
allData_dt <- get(load(file.path(inFolder, hicds, paste0(exprds, "_", familyData), familySubData, "allData_dt.Rdata")))

allData_dt <- allData_dt[allData_dt$dist <= maxDist,]

allData_dt$sameTAD <- ifelse(allData_dt$sameTAD == 1, "sameTAD",
                             ifelse(allData_dt$sameTAD == 0, "diffTAD", NA))
stopifnot(!is.na(allData_dt$sameTAD))

#allData_dt$sameFamily <- ifelse(allData_dt$sameFamily == 1, "sameFam",
#                             ifelse(allData_dt$sameFamily == 0, "diffFam", NA))
#stopifnot(!is.na(allData_dt$sameFamily))

dist_histBreaks_vect <- seq(0, maxDist, length.out=10+1)
dist_histBreaks_vect_labs <- paste0("]", dist_histBreaks_vect[1:(length(dist_histBreaks_vect)-1)]/1000, ", ",dist_histBreaks_vect[2:length(dist_histBreaks_vect)]/1000, "]")
dist_histBreaks_vect_labs[1] <- sub("]", "[", dist_histBreaks_vect_labs[1])

allData_dt$dist_cat <- foreach(i = 1:nrow(allData_dt), .combine='c' ) %dopar% { which(hist(allData_dt$dist[i], breaks=dist_histBreaks_vect, plot=FALSE)$counts == 1) }

outFile <- file.path(outFolder, "allData_dt.Rdata")
save(allData_dt, file=outFile, version=2)  
cat(paste0("... written: ", outFile, "\n"))

#allData_dt$cond <- interaction(allData_dt$sameTAD, allData_dt$sameFamily)
allData_dt$cond <- interaction(allData_dt$sameTAD)

#allData_dt$cond <- factor(allData_dt$cond, levels = c( "sameTAD.sameFam", "sameTAD.diffFam",   "diffTAD.sameFam", "diffTAD.diffFam"))
allData_dt$cond <- factor(allData_dt$cond, levels = c( "sameTAD",   "diffTAD"))
stopifnot(!is.na(allData_dt$cond))

allData_dt$dist_cat <- factor(allData_dt$dist_cat, levels=as.character(sort(unique(allData_dt$dist_cat))))

subTit <- paste0("max dist. = ", maxDist/1000, " kb")
my_ylab <- "Gene1-gene2 expr. corr."
my_xlab <- " Gene1-gene2 dist. range [kb]"

p_coexpr_boxplot <- ggplot(allData_dt, aes(x=dist_cat, y = coexpr, fill = cond, color = cond)) + 
 
  geom_boxplot(notch = TRUE, outlier.shape=NA)+
  ggtitle(paste0(hicds, " - ", exprds), subtitle = paste0(subTit))+
  scale_x_discrete(name=my_xlab, labels =dist_histBreaks_vect_labs )+
  scale_y_continuous(name=paste0(my_ylab),
                     breaks = scales::pretty_breaks(n = 20))+
  
  labs(fill  = paste0(""),  color=paste0("")) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    axis.line.x= element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
    axis.text.x =element_text(color="black", hjust=0.5,vjust = 0.5, size=10, face="bold"),
    # axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.text = element_text(size=12),
    legend.key = element_blank(),
    
    legend.key.size = unit(1.2, 'cm'),
    
    legend.title = element_text(face="bold", size=12)
  )
  
outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_sameTAD_diffTAD_boxplot.", plotType))
ggsave(plot = p_coexpr_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


p_coexpr_boxplot_jitter <- p_coexpr_boxplot +
  geom_point( position=position_jitterdodge(), stroke=0.8, shape=21, alpha=0.8) 

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_sameTAD_diffTAD_boxplot_jitter.", plotType))
ggsave(plot = p_coexpr_boxplot_jitter, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_sameTAD_diffTAD_countPairs_dt.txt"))
write.table(table(allData_dt$cond, allData_dt$dist_cat), file =outFile, col.names=T, row.names = T, quote=F, sep="\t")
cat(paste0("... written: ", outFile, "\n"))

count_dt <- as.data.frame(table(allData_dt$cond, allData_dt$dist_cat))
# count_dt$Freq_log10 <- log10(count_dt$Freq)

my_ylab <- "# gene pairs"

p_count_barplot <- ggplot(count_dt, aes(x=Var2, y = Freq, fill = Var1, color = Var1)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle(paste0(hicds, " - ", exprds), subtitle = paste0(subTit))+
  scale_x_discrete(name=my_xlab, labels =dist_histBreaks_vect_labs )+
  scale_y_continuous(name=paste0(my_ylab),
                     breaks = scales::pretty_breaks(n = 20))+
  
  labs(fill  = paste0(""),  color=paste0("")) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    axis.line.x= element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
    axis.text.x =element_text(color="black", hjust=0.5,vjust = 0.5, size=10, face="bold"),
    # axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.text = element_text(size=12),
    legend.key = element_blank(),
    
    legend.key.size = unit(1, 'cm'),
    
    legend.title = element_text(face="bold", size=12)
  )

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_sameTAD_diffTAD_countPairs_barplot.", plotType))
ggsave(plot = p_count_barplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))



######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





# foo_breaks <- c(0,10,20,30,40,50)
# hist(5, breaks=foo_breaks)$breaks
# # c(0,10,20,30,40,50)
# hist(5, breaks=foo_breaks)$counts
# # 1 0 0 0 0
# which(hist(5, breaks=foo_breaks)$counts == 1)
# # 1
# which(hist(0, breaks=foo_breaks)$counts == 1)
# # 1
# which(hist(10, breaks=foo_breaks)$counts == 1)
# # 1
# which(hist(50, breaks=foo_breaks)$counts == 1)
# # 5
# which(hist(39, breaks=foo_breaks)$counts == 1)
# # 4
# which(hist(40, breaks=foo_breaks)$counts == 1)
# # 4
# which(hist(41, breaks=foo_breaks)$counts == 1)
# # 5
# which(hist(11, breaks=foo_breaks)$counts == 1)
# # 2
# which(hist(55, breaks=foo_breaks)$counts == 1)
# # ERROR
# which(hist(-2, breaks=foo_breaks)$counts == 1)
# # ERROR

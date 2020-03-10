startTime <- Sys.time()
cat(paste0("> Rscript coexprDist_sameTAD_boxplot_cmp_permut_top50.R\n"))

# Rscript coexprDist_sameTAD_boxplot_cmp_permut_top50.R

options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

registerDoMC(40)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

hicds = "ENCSR489OCU_NCI-H460_40kb"
exprds= "TCGAluad_norm_luad"

nTop <- 50

args <- commandArgs(trailingOnly = TRUE)

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  # all_hicds <- all_hicds[grepl("_RANDOMMIDPOS_", all_hicds) ]
  all_hicds <- all_hicds[! (grepl("_RANDOM", all_hicds) | grepl("_PERMUT", all_hicds)) ]
  # all_hicds <- all_hicds[!grepl("_RANDOMMIDPOSDISC", all_hicds)]
  # all_hicds <- "ENCSR489OCU_NCI-H460_40kb"
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  stopifnot(length(args) == 2)
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}

plotType <- "svg"
myHeightGG <- myWidthGG <- 7

plotType <- "png"
myHeightGG <- myWidthGG <- 7


familyData <- "hgnc"
familySubData <- "hgnc_family_short"

inFolder <- file.path("AUC_COEXPRDIST_WITHFAM_SORTNODUP")

maxDist <- 500*1000

rd_patterns <- paste0("", c("PERMUTG2T", "RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMSHIFT" , "RANDOMMIDPOSDISC"))
rd_patterns <- c( "RANDOMMIDPOS",  "RANDOMMIDPOSDISC")
rd_patt <- rd_patterns[1]

dist_histBreaks_vect <- seq(0, maxDist, length.out=10+1)
dist_histBreaks_vect_labs <- paste0("]", dist_histBreaks_vect[1:(length(dist_histBreaks_vect)-1)]/1000, ", ",dist_histBreaks_vect[2:length(dist_histBreaks_vect)]/1000, "]")
dist_histBreaks_vect_labs[1] <- sub("]", "[", dist_histBreaks_vect_labs[1])

tad_signifThresh <- 0.01

hicds = all_hicds[1]

all_signif_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    outFolder <- file.path(paste0("COEXPRDIST_SAMETAD_BOXPLOT_CMP_PERMUT_TOP", nTop), hicds, exprds)
    dir.create(outFolder, recursive = TRUE)
    
    combPval_file <- file.path(pipOutFolder, hicds, exprds, "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")
    stopifnot(file.exists(combPval_file))
    pval <- get(load(combPval_file))
    adjPval <- p.adjust(pval, method="BH")
    adjPval <- sort(adjPval)
    top_tads <- names(adjPval)[1:nTop]
    signif_tads <- names(adjPval)[adjPval <= tad_signifThresh]
    
    
    g2t_dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t <- setNames(g2t_dt$region, g2t_dt$entrezID)

        
    cat(paste0("... start : ", hicds, " - ", exprds, "\n"))
    
    
    allData_dt <- get(load(file.path(inFolder, hicds, paste0(exprds, "_", familyData), familySubData, "allData_dt.Rdata")))
    
    allData_dt <- allData_dt[allData_dt$dist <= maxDist,]
    
    allData_dt$sameTAD <- ifelse(allData_dt$sameTAD == 1, "sameTAD",
                                 ifelse(allData_dt$sameTAD == 0, "diffTAD", NA))
    stopifnot(!is.na(allData_dt$sameTAD))
    
    allData_dt$sameFamily <- NULL

    
    to_merge_dt <- allData_dt[ allData_dt$sameTAD == "sameTAD",]
    stopifnot(nrow(to_merge_dt) > 0)
    
    to_merge_dt$gene1_tad <- g2t[as.character(to_merge_dt$gene1)]
    stopifnot(!is.na(to_merge_dt$gene1_tad))
    
    to_merge_dt$gene2_tad <- g2t[as.character(to_merge_dt$gene2)]
    stopifnot(!is.na(to_merge_dt$gene2_tad))
    
    stopifnot(to_merge_dt$gene1_tad == to_merge_dt$gene2_tad)
    
    to_merge_dt <- to_merge_dt[to_merge_dt$gene1_tad %in% top_tads,]
    
    to_merge_dt$dist_cat <- foreach(i = 1:nrow(to_merge_dt), .combine='c' ) %dopar% { which(hist(to_merge_dt$dist[i], breaks=dist_histBreaks_vect, plot=FALSE)$counts == 1) }
    
    to_merge_dt$dataType <- "OBSERVED"
    
    signif_to_merge_dt <- to_merge_dt[to_merge_dt$gene1_tad %in% signif_tads,]
    
    signif_merged_dt <- signif_to_merge_dt
    merged_dt <- to_merge_dt
    
    rm("top_tads")
    rm("signif_tads")
    rm("g2t")
    rm("adjPval")
    
    for(rd_patt in rd_patterns) {
      
      cat(paste0("... start ", rd_patt, "\n"))
      
      rd_hicds <- gsub("_40kb", paste0("_", rd_patt, "_40kb"), hicds)
      
      rd_combPval_file <- file.path(pipOutFolder, rd_hicds, exprds, "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")
      stopifnot(file.exists(rd_combPval_file))
      rd_pval <- get(load(rd_combPval_file))
      rd_adjPval <- p.adjust(rd_pval, method="BH")
      rd_adjPval <- sort(rd_adjPval)
      rd_top_tads <- names(rd_adjPval)[1:nTop]
      rd_signif_tads <- names(rd_adjPval)[rd_adjPval <= tad_signifThresh]
      
      rd_g2t_dt <- read.delim(file.path(rd_hicds, "genes2tad", "all_genes_positions.txt"), header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
      rd_g2t <- setNames(rd_g2t_dt$region, rd_g2t_dt$entrezID)
      
      
      allData_dt <- get(load(file.path(inFolder, rd_hicds, paste0(exprds, "_", familyData), familySubData, "allData_dt.Rdata")))
      
      allData_dt <- allData_dt[allData_dt$dist <= maxDist,]
      
      allData_dt$sameTAD <- ifelse(allData_dt$sameTAD == 1, "sameTAD",
                                   ifelse(allData_dt$sameTAD == 0, "diffTAD", NA))
      stopifnot(!is.na(allData_dt$sameTAD))
      
      allData_dt$sameFamily <- NULL
      
      
      to_merge_dt <- allData_dt[allData_dt$sameTAD == "sameTAD",]
      stopifnot(nrow(to_merge_dt) > 0)
      
      to_merge_dt$gene1_tad <- rd_g2t[as.character(to_merge_dt$gene1)]
      stopifnot(!is.na(to_merge_dt$gene1_tad))
      
      to_merge_dt$gene2_tad <- rd_g2t[as.character(to_merge_dt$gene2)]
      stopifnot(!is.na(to_merge_dt$gene2_tad))
      
      stopifnot(to_merge_dt$gene1_tad == to_merge_dt$gene2_tad)
      
      to_merge_dt <- to_merge_dt[to_merge_dt$gene1_tad %in% rd_top_tads,]
      
      to_merge_dt$dist_cat <- foreach(i = 1:nrow(to_merge_dt), .combine='c' ) %dopar% { which(hist(to_merge_dt$dist[i], breaks=dist_histBreaks_vect, plot=FALSE)$counts == 1) }
      
      to_merge_dt$dataType <- rd_patt
      
      merged_dt <- rbind(merged_dt, to_merge_dt)
      
      signif_to_merge_dt <- to_merge_dt[to_merge_dt$gene1_tad %in% rd_signif_tads,]
      signif_merged_dt <- rbind(signif_merged_dt, signif_to_merge_dt)
      
    }
    
    rm("allData_dt")
    
    outFile <- file.path(outFolder, "merged_dt.Rdata")
    save(merged_dt, file=outFile, version=2)  
    cat(paste0("... written: ", outFile, "\n"))


    stopifnot(merged_dt$sameTAD == "sameTAD")
    
    merged_dt$dataType <- factor(merged_dt$dataType, levels = c( "OBSERVED", rd_patterns))
    stopifnot(!is.na(merged_dt$dataType))
    
    merged_dt$dist_cat <- factor(merged_dt$dist_cat, levels=as.character(sort(unique(merged_dt$dist_cat))))
    
    subTit <- paste0("sameTAD gene pairs; max dist. = ", maxDist/1000, " kb; top ", nTop, " TADs")
    my_ylab <- "Gene1-gene2 expr. corr."
    my_xlab <- " Gene1-gene2 dist. range [kb]"
    
    p_coexpr_boxplot <- ggplot(merged_dt, aes(x=dist_cat, y = coexpr, fill = dataType, color = dataType)) + 
      
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
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_cmp_obs_permut_sameTAD_boxplot.", plotType))
    ggsave(plot = p_coexpr_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    p_coexpr_boxplot_jitter <- p_coexpr_boxplot +
      geom_point( position=position_jitterdodge(), stroke=0.8, shape=21, alpha=0.8) 
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_cmp_obs_permut_sameTAD_boxplot_jitter.", plotType))
    ggsave(plot = p_coexpr_boxplot_jitter, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_cmp_obs_permut_sameTAD_countPairs_dt.txt"))
    write.table(table(merged_dt$dataType, merged_dt$dist_cat), file =outFile, col.names=T, row.names = T, quote=F, sep="\t")
    cat(paste0("... written: ", outFile, "\n"))
    
    count_dt <- as.data.frame(table(merged_dt$dataType, merged_dt$dist_cat))
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
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_cmp_obs_permut_sameTAD_countPairs_barplot.", plotType))
    ggsave(plot = p_count_barplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
    cat(paste0("... written: ", outFile, "\n"))
    
    signif_merged_dt$hicds <- hicds
    signif_merged_dt$exprds <- exprds
    signif_merged_dt
  } # end foreach exprds
  exprds_dt
} # end foreach hicds
   
save(all_signif_dt, file="all_signif_dt.Rdata", version=2)

stopifnot(all_signif_dt$sameTAD == "sameTAD")

all_signif_dt$dataType <- factor(all_signif_dt$dataType, levels = c( "OBSERVED", rd_patterns))
stopifnot(!is.na(all_signif_dt$dataType))

all_signif_dt$dist_cat <- factor(all_signif_dt$dist_cat, levels=as.character(sort(unique(all_signif_dt$dist_cat))))

subTit <- paste0("sameTAD gene pairs; max dist. = ", maxDist/1000, " kb; top ", nTop, " TADs")
my_ylab <- "Gene1-gene2 expr. corr."
my_xlab <- " Gene1-gene2 dist. range [kb]"


outFolder <- file.path(paste0("COEXPRDIST_SAMETAD_BOXPLOT_CMP_PERMUT_TOP", nTop))

nDS <- length(unique(file.path(all_signif_dt$hicds, all_signif_dt$exprds)))

subTit <- paste0("sameTAD gene pairs; max dist. = ", maxDist/1000, " kb; signif. TADs (adj. p-val <=", tad_signifThresh, ")")
plotTit <-  paste0("all DS (n=", nDS, ")")
my_ylab <- "Gene1-gene2 expr. corr."
my_xlab <- " Gene1-gene2 dist. range [kb]"

all_coexpr_boxplot <- ggplot(all_signif_dt, aes(x=dist_cat, y = coexpr, fill = dataType, color = dataType)) + 
  
  geom_boxplot(notch = TRUE, outlier.shape=NA)+
  ggtitle(paste0(plotTit), subtitle = paste0(subTit))+
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

outFile <- file.path(outFolder, paste0("allDS_coexpr_cmp_obs_permut_sameTAD_signifTADs_boxplot.", plotType))
ggsave(plot = all_coexpr_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


count_dt <- as.data.frame(table(all_signif_dt$dataType, all_signif_dt$dist_cat))
# count_dt$Freq_log10 <- log10(count_dt$Freq)

my_ylab <- "# gene pairs"

all_count_barplot <- ggplot(count_dt, aes(x=Var2, y = Freq, fill = Var1, color = Var1)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle(plotTit, subtitle = paste0(subTit))+
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

outFile <- file.path(outFolder, paste0("allDS_coexpr_cmp_obs_permut_sameTAD_signifTADs_countPairs_barplot.", plotType))
ggsave(plot = all_count_barplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
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

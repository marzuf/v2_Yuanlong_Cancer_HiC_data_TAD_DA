startTime <- Sys.time()
cat(paste0("> Rscript coexprDist_sameFam_sameTAD_boxplot_cmp_permut.R\n"))

# Rscript coexprDist_sameFam_sameTAD_boxplot_cmp_permut.R
# Rscript coexprDist_sameFam_sameTAD_boxplot_cmp_permut.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad

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
rd_patterns <- c( "RANDOMMIDPOS",  "RANDOMMIDPOSDISC", "RANDOMMIDPOSSTRICT")
rd_patt <- rd_patterns[1]

dist_histBreaks_vect <- seq(0, maxDist, length.out=10+1)
dist_histBreaks_vect_labs <- paste0("]", dist_histBreaks_vect[1:(length(dist_histBreaks_vect)-1)]/1000, ", ",dist_histBreaks_vect[2:length(dist_histBreaks_vect)]/1000, "]")
dist_histBreaks_vect_labs[1] <- sub("]", "[", dist_histBreaks_vect_labs[1])

hicds = all_hicds[1]

FOO <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  FOO <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    cat(paste0("... start : ", hicds, " - ", exprds, "\n"))
    
    outFolder <- file.path("COEXPRDIST_SAMEFAM_SAMETAD_BOXPLOT_CMP_PERMUT", hicds, exprds)
    dir.create(outFolder, recursive = TRUE)
    
    outFile_model <- file.path(outFolder, "allData_lm_summary.txt")
    file.remove(outFile_model)
    
    
    allData_dt <- get(load(file.path(inFolder, hicds, paste0(exprds, "_", familyData), familySubData, "allData_dt.Rdata")))
    
    allData_dt <- allData_dt[allData_dt$dist <= maxDist,]
    
    allData_dt$sameTAD <- ifelse(allData_dt$sameTAD == 1, "sameTAD",
                                 ifelse(allData_dt$sameTAD == 0, "diffTAD", NA))
    stopifnot(!is.na(allData_dt$sameTAD))
    
    allData_dt$sameFamily <- ifelse(allData_dt$sameFamily == 1, "sameFam",
                                    ifelse(allData_dt$sameFamily == 0, "diffFam", NA))
    stopifnot(!is.na(allData_dt$sameFamily))
    
    m1 <- lm(coexpr ~ sameTAD, data = allData_dt)
    m2 <- lm(coexpr ~ sameFamily, data = allData_dt)
    m3 <- lm(coexpr ~ sameTAD+sameFamily, data = allData_dt)
    
    sink(outFile_model, append=TRUE)
    print("************************** OBSERVED\n")
    print(summary(m1))
    print("\n----------------------------------------------------\n")
    print(summary(m2))
    print("\n----------------------------------------------------\n")
    print(summary(m3))
    sink()
    
    to_merge_dt <- allData_dt[allData_dt$sameFamily == "sameFam" & allData_dt$sameTAD == "sameTAD",]
    stopifnot(nrow(to_merge_dt) > 0)
    to_merge_dt$dist_cat <- foreach(i = 1:nrow(to_merge_dt), .combine='c' ) %dopar% { which(hist(to_merge_dt$dist[i], breaks=dist_histBreaks_vect, plot=FALSE)$counts == 1) }
    
    to_merge_dt$dataType <- "OBSERVED"
    
    
    merged_dt <- to_merge_dt
    
    
    for(rd_patt in rd_patterns) {
      
      cat(paste0("... start ", rd_patt, "\n"))
      
      rd_hicds <- gsub("_40kb", paste0("_", rd_patt, "_40kb"), hicds)
      
      allData_dt <- get(load(file.path(inFolder, rd_hicds, paste0(exprds, "_", familyData), familySubData, "allData_dt.Rdata")))
      
      allData_dt <- allData_dt[allData_dt$dist <= maxDist,]
      
      allData_dt$sameTAD <- ifelse(allData_dt$sameTAD == 1, "sameTAD",
                                   ifelse(allData_dt$sameTAD == 0, "diffTAD", NA))
      stopifnot(!is.na(allData_dt$sameTAD))
      
      allData_dt$sameFamily <- ifelse(allData_dt$sameFamily == 1, "sameFam",
                                      ifelse(allData_dt$sameFamily == 0, "diffFam", NA))
      stopifnot(!is.na(allData_dt$sameFamily))
      
      
      m1 <- lm(coexpr ~ sameTAD, data = allData_dt)
      m2 <- lm(coexpr ~ sameFamily, data = allData_dt)
      m3 <- lm(coexpr ~ sameTAD+sameFamily, data = allData_dt)
      
      sink(outFile_model, append=TRUE)
      print(paste0("************************** ", rd_patt, "\n"))
      print(summary(m1))
      print("\n----------------------------------------------------\n")
      print(summary(m2))
      print("\n----------------------------------------------------\n")
      print(summary(m3))
      sink()
      
      to_merge_dt <- allData_dt[allData_dt$sameFamily == "sameFam" & allData_dt$sameTAD == "sameTAD",]
      stopifnot(nrow(to_merge_dt) > 0)
      to_merge_dt$dist_cat <- foreach(i = 1:nrow(to_merge_dt), .combine='c' ) %dopar% { which(hist(to_merge_dt$dist[i], breaks=dist_histBreaks_vect, plot=FALSE)$counts == 1) }
      
      to_merge_dt$dataType <- rd_patt
      
      merged_dt <- rbind(merged_dt, to_merge_dt)
      
    }
    
    rm("allData_dt")
    
    outFile <- file.path(outFolder, "merged_dt.Rdata")
    save(merged_dt, file=outFile, version=2)  
    cat(paste0("... written: ", outFile, "\n"))
    
    stopifnot(merged_dt$sameFamily == "sameFam")
    stopifnot(merged_dt$sameTAD == "sameTAD")
    
    merged_dt$dataType <- factor(merged_dt$dataType, levels = c( "OBSERVED", rd_patterns))
    stopifnot(!is.na(merged_dt$dataType))
    
    merged_dt$dist_cat <- factor(merged_dt$dist_cat, levels=as.character(sort(unique(merged_dt$dist_cat))))
    
    subTit <- paste0("sameFam+sameTAD gene pairs; max dist. = ", maxDist/1000, " kb")
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
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_cmp_obs_permut_sameFam_sameTAD_boxplot.", plotType))
    ggsave(plot = p_coexpr_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    p_coexpr_boxplot_jitter <- p_coexpr_boxplot +
      geom_point( position=position_jitterdodge(), stroke=0.8, shape=21, alpha=0.8) 
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_cmp_obs_permut_sameFam_sameTAD_boxplot_jitter.", plotType))
    ggsave(plot = p_coexpr_boxplot_jitter, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_cmp_obs_permut_sameFam_sameTAD_countPairs_dt.txt"))
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
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_cmp_obs_permut_sameFam_sameTAD_countPairs_barplot.", plotType))
    ggsave(plot = p_count_barplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    cat(paste0("... written: ", outFile_model, "\n"))
    
    
    
    

    
    
    
  } # end foreach exprds
} # end foreach hicds
   



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

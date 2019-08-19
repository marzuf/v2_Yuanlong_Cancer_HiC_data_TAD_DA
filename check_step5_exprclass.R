
options(scipen=100)

# Rscript check_step5_exprclass.R

script_name <- "check_step5_exprclass.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)

registerDoMC(10)

# mainFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
mainFolder <- "."
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")

hicds <- "GSE105318_DLD1_40kb"
exprds <- "TCGAcoad_msi_mss"

 

script0_name <- "0_prepGeneData"
script5_name <- "5_runPermutationsMedian"

require(ggplot2)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- "CHECK_STEP5_EXPRCLASS"
dir.create(outFolder, recursive = TRUE)

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myWidthGG <- 12
myHeightGG <- 8

maxPermutToTake <- 1000


all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds


all_data_list <- foreach(hicds = all_hicds) %dopar% {
  

  g2t_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
  stopifnot(file.exists(g2t_file))
  g2t_dt <- read.delim(g2t_file, header=F, col.names=c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
  g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
  
  
  exprds = all_exprds[[paste0(hicds)]][1]
  exprds_list <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
    
    
    permDT_file <- file.path(pipFolder, hicds, exprds, script5_name, "permutationsDT.Rdata")
    stopifnot(file.exists(permDT_file))
    permDT <- get(load(permDT_file))
    nPermut <- min(ncol(permDT), maxPermutToTake)
    permDT <- permDT[,1:nPermut]
    
    fpkm_file <- file.path(pipFolder, hicds, exprds, script0_name, "rna_fpkmDT.Rdata")
    stopifnot(file.exists(fpkm_file))
    fpkmDT <- get(load(fpkm_file))
    medianExpr <- apply(fpkmDT, 1, median)
    stopifnot(!is.na(medianExpr))
    fpkmDT$medianExpr <- medianExpr
    
    geneList_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneList_file))
    geneList <- get(load(geneList_file))
    
    stopifnot(setequal(rownames(permDT), geneList))
    
    stopifnot(names(geneList) %in% rownames(fpkmDT))
    
    fpkmDT <- fpkmDT[rownames(fpkmDT) %in% names(geneList),]
    
    fpkmDT$entrezID <- sapply(rownames(fpkmDT), function(x) as.character(geneList[names(geneList) == x]))
    stopifnot(!is.na(fpkmDT$entrezID))
    stopifnot(is.character(fpkmDT$entrezID))
    stopifnot(fpkmDT$entrezID %in% g2t_dt$entrezID)
    fpkmDT <- fpkmDT[,c("medianExpr", "entrezID")]
    rownames(fpkmDT) <- NULL
    obs_fpkmDT <- merge(fpkmDT, g2t_dt[, c("entrezID", "region")], all.x=TRUE, all.y=FALSE, by="entrezID")
    stopifnot(nrow(obs_fpkmDT) == nrow(fpkmDT))
    stopifnot(!is.na(obs_fpkmDT$region))
    
    obs_TAD_meanExpr_DT <- aggregate(medianExpr ~ region, data=obs_fpkmDT, FUN=mean)
    stopifnot(!is.na(obs_TAD_meanExpr_DT$medianExpr))
    obs_TAD_meanExpr_DT <- obs_TAD_meanExpr_DT[order(obs_TAD_meanExpr_DT$medianExpr),]
    obs_TAD_meanExpr_DT$region <- factor(as.character(obs_TAD_meanExpr_DT$region), levels = as.character(obs_TAD_meanExpr_DT$region)) 
    
    all_regs <- as.character(obs_TAD_meanExpr_DT$region)
    
    # do the same for each permut
    meanExpr_permDT <- foreach(i_perm = 1:ncol(permDT), .combine='rbind') %dopar% {
      
      assignDT <- data.frame(entrezID = rownames(permDT), region = permDT[,i_perm])
      
      perm_exprDT <- merge(assignDT, fpkmDT, by="entrezID",all.x=TRUE,all.y=FALSE)
      stopifnot(nrow(perm_exprDT) == nrow(assignDT))
      stopifnot(!is.na(perm_exprDT$medianExpr))
      
      perm_TAD_mean <- aggregate(medianExpr ~ region, data=perm_exprDT, FUN=mean)
      stopifnot(setequal(all_regs, perm_TAD_mean$region))
      
      perm_TAD_mean
    }
    meanExpr_permDT$region <- factor(as.character(meanExpr_permDT$region), levels = all_regs)
    colnames(meanExpr_permDT)[colnames(meanExpr_permDT) == "medianExpr"] <- "TAD_mean_medianExpr"
    
    TAD_mean_medianExpr_log10 <- log10(meanExpr_permDT$TAD_mean_medianExpr)
    
    p_var <-  ggplot(meanExpr_permDT, aes_string(x = "region", y = "TAD_mean_medianExpr_log10")) + 
      geom_boxplot()+
      # coord_cartesian(expand = FALSE) +
      ggtitle(paste0("TAD_mean_medianExpr "), subtitle = paste0(nPermut, " permut."))+
      scale_x_discrete(name="TAD ranked by obs. mean medianExpr.")+
      scale_y_continuous(name=paste0("permut. TAD_mean_medianExpr"),
                         breaks = scales::pretty_breaks(n = 20))+
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
        # axis.text.x =  element_text(color="black", hjust=1,vjust = 0.5, angle=90),
        axis.text.x =  element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        legend.background =  element_rect(),
        legend.key = element_blank(),
        legend.title = element_text(face="bold")
      )
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_TAD", "_", "mean_medianExpr", "_", nPermut, "permut_boxplot.", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    meanPermut <- aggregate(TAD_mean_medianExpr ~ region, FUN=mean, data=meanExpr_permDT)
    stopifnot(setequal(meanPermut$region, obs_TAD_meanExpr_DT$region))
    stopifnot(nrow(meanPermut) == nrow(obs_TAD_meanExpr_DT))
    permutValues <- setNames(meanPermut$TAD_mean_medianExpr, meanPermut$region)
    obsValues <- setNames(obs_TAD_meanExpr_DT$medianExpr, obs_TAD_meanExpr_DT$region)
    stopifnot(setequal(names(permutValues), all_regs))
    stopifnot(setequal(names(obsValues), all_regs))
    
    myx <- log10(obsValues[all_regs])
    myy <- log10(permutValues[all_regs])
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_TAD", "_", "mean_medianExpr", "_mean", nPermut, "permut_densplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x=myx,
      y=myy,
      xlab=paste0("obs. values (log10)"),
      ylab=paste0("mean permut. values (log10)"),
      main = "TAD mean median expression"
    )
    addCorr(x=myx, y = myy, legPos = "topleft", bty="n")
    mtext(text = paste0("nPermut=", nPermut), side=3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))

  } # end-foreach iterating over exprds
} # end-foreach iterating over hicds
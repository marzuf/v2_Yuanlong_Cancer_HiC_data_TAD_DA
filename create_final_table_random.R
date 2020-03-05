### CREATE TABLE
# -----------------------
#  TAD_ID | genes | meanFC | meanCorr | adj. comb. p-value | signif. at FDR 0.1 ?  | signif. at FDR 0.2  

# Rscript create_final_table_random.R

# 12.08.2019 -> added start and end positions

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "create_final_table_random.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script11same_name <- "11sameNbr_runEmpPvalCombined"
script19_name <- "19onlyFC_SAM_emp_measurement" # 
script19sameNbr_name <- "19sameNbr_SAM_emp_measurement"
script8c_name <- "8cOnlyRatioDownFastSave_runAllDown" # -> 

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
require(reshape2)

pval001_col <- "dodgerblue3"
pval005_col <- "darkorange3"
fdr01_col <- "goldenrod"
fdr02_col <- "darkolivegreen"
fdr_pval_col <- "indianred4"


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myWidthGG <- 12
myHeightGG <- 8


pipFolder<- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "CREATE_FINAL_TABLE_RANDOM" # "CREATE_FINAL_TABLE_10000" for 10000 permut
dir.create(outFolder, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
hicds <- args[1]
exprds <- args[2]

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
entrez2symb_dt$symbol <- as.character(entrez2symb_dt$symbol)
stopifnot(!duplicated(entrez2symb_dt$entrezID))


fdr_thresh1 <- 0.1
fdr_thresh2 <- 0.2

buildTable <- TRUE


if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  all_hicds <- all_hicds[grepl("_RANDOMMIDPOS_", all_hicds) ]
  # all_hicds <- all_hicds[grepl("_RANDOM", all_hicds) | grepl("_PERMUT", all_hicds) ]
  # all_hicds <- all_hicds[!grepl("_RANDOMMIDPOSDISC", all_hicds)]
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}

minTADsize <- 3


hicds = all_hicds[1]

if(buildTable) {
  
  
  all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %do% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start building DT for :", hicds, " - ", exprds, "\n")
      
      
      ### RETRIEVE THE GENE2TAD ASSIGNMENT
      g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      

      ### RETRIEVE THE TAD ASSIGNMENT
      tadFile <- file.path(pipFolder, hicds, "genes2tad", "all_assigned_regions.txt")
      stopifnot(file.exists(tadFile))
      tad_DT <- read.delim(tadFile, header=F, col.names = c("chromo","region", "start", "end"), stringsAsFactors = FALSE)
      tad_DT <- tad_DT[grepl("_TAD", tad_DT$region),]
      stopifnot(nrow(tad_DT) > 0)

      ### RETRIEVE THE GENES USED IN THE PIPELINE - script0
      script0_name <- "0_prepGeneData"
      stopifnot(dir.exists(file.path(pipOutFolder, hicds, exprds, script0_name)))
      geneListFile <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneListFile))
      pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted
      stopifnot(pipeline_geneList %in% g2t_DT$entrezID)    
      g2t_DT <- g2t_DT[g2t_DT$entrezID %in% pipeline_geneList,]


      # RETRIEVE FDR DATA FOR MEAN CORR
      # the same for meanCorr
      meanCorr_FDR_file <-  file.path(pipOutFolder, hicds, exprds, script19sameNbr_name, "meanCorr_empFDR.Rdata")
      stopifnot(file.exists(meanCorr_FDR_file))
      all_corr_FDR <- eval(parse(text = load(meanCorr_FDR_file)))

      fcEmpFDR_file <- file.path(pipOutFolder, hicds, exprds, script19_name, "empFDR_list.Rdata")
      stopifnot(file.exists(fcEmpFDR_file))
      fc_empFDR <- eval(parse(text = load(fcEmpFDR_file)))
      
      fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
      stopifnot(file.exists(fc_file))
      tad_fc <- eval(parse(text = load(fc_file)))
      all_regs <- names(tad_fc)

      stopifnot(all_regs %in% g2t_DT$region)
      stopifnot(all_regs %in% tad_DT$region)
      all_regs_start <- setNames(tad_DT$start, tad_DT$region)
      all_regs_end <- setNames(tad_DT$end, tad_DT$region)

      
      corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
      stopifnot(file.exists(corr_file))  
      tad_corr <- eval(parse(text = load(corr_file)))
      stopifnot(setequal(all_regs, names(tad_corr)))

      comb_empPval_file <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
      stopifnot(file.exists(comb_empPval_file))
      comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
      stopifnot(setequal(all_regs, names(comb_empPval)))
      comb_empPval <- comb_empPval[all_regs]
      # ADJUST THE PVAL
      tad_adjCombPval <- p.adjust(comb_empPval, method="BH")
      stopifnot(names(tad_adjCombPval) == all_regs)


      fc_FDR <- fc_empFDR[["empFDR_logFC"]]
      # the 1st FC value for which FDR is smaller than the thresh
      fc_FDR_co1 <- min(as.numeric(as.character(na.omit(names(fc_FDR)[fc_FDR <= fdr_thresh1]))))
      # will return Inf is none below the thresh
      fc_FDR_co2 <- min(as.numeric(as.character(na.omit(names(fc_FDR)[fc_FDR <= fdr_thresh2]))))
      
      corr_FDR <- all_corr_FDR[["empFDR"]]
      corr_FDR_co1 <- min(as.numeric(as.character(na.omit(names(corr_FDR)[corr_FDR <= fdr_thresh1])))) # will return Inf is none below the thresh
      corr_FDR_co2 <- min(as.numeric(as.character(na.omit(names(corr_FDR)[corr_FDR <= fdr_thresh2])))) # will return Inf is none below the thresh
      
      FDR_signif_co1 <- setNames(abs(tad_fc[all_regs]) >= fc_FDR_co1 &
                                   tad_corr[all_regs] >= corr_FDR_co1, all_regs)
      
      stopifnot(setequal(all_regs, names(FDR_signif_co1)))
      
      
      FDR_signif_co2 <- setNames(abs(tad_fc[all_regs]) >= fc_FDR_co2 &
                                   tad_corr[all_regs] >= corr_FDR_co2, all_regs)
      
      
      
      stopifnot(setequal(all_regs, names(FDR_signif_co2)))
      
      
      # RETRIEVE RATIO DOWN
      rd_file <- file.path(pipOutFolder, hicds, exprds, script8c_name, "all_obs_ratioDown.Rdata")
      stopifnot(file.exists(rd_file))
      tad_rD <- eval(parse(text = load(rd_file)))
      stopifnot(setequal(names(tad_rD), all_regs))
      
      
      tad_genes <- sapply(all_regs, function(x) {
        x_genes <- g2t_DT$entrezID[g2t_DT$region == x]
        stopifnot(length(x_genes) >= minTADsize)
        stopifnot(x_genes %in% entrez2symb_dt$entrezID)
        x_symbols <- entrez2symb_dt$symbol[entrez2symb_dt$entrezID %in% x_genes]
        stopifnot(length(x_symbols) == length(x_genes))
        x_symbols <- sort(x_symbols)
        paste0(x_symbols, collapse=",")
      })
      names(tad_genes) <- all_regs

      stopifnot(all_regs %in% names(all_regs_end))
      stopifnot(all_regs %in% names(all_regs_start))
      
      
      exprds_hicds_dt <- data.frame(
        hicds = hicds,
        exprds = exprds,
        region = all_regs,

        start = all_regs_start[all_regs],
        end = all_regs_end[all_regs],

        region_genes = tad_genes[all_regs],
        meanLogFC = tad_fc[all_regs],
        meanCorr = tad_corr[all_regs],
        
        ratioDown = tad_rD[all_regs],
        
        meanLogFC_FDR_co1 = fc_FDR_co1,
        meanLogFC_FDR_co2 = fc_FDR_co2,
        
        meanCorr_FDR_co1 = corr_FDR_co1,
        meanCorr_FDR_co2 = corr_FDR_co2,
        
        
        adjPvalComb = tad_adjCombPval[all_regs],
        signifFDRco1 = FDR_signif_co1[all_regs],
        signifFDRco2 = FDR_signif_co2[all_regs],
        stringsAsFactors = FALSE
      )
      rownames(exprds_hicds_dt) <- NULL
      colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "signifFDRco1"] <- paste0("signifFDR_",
                                                                                       fdr_thresh1)
      colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "signifFDRco2"] <- paste0("signifFDR_",
                                                                                       fdr_thresh2)
      
      colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "meanLogFC_FDR_co1"] <- paste0("meanLogFC_thresh_FDR",
                                                                                            fdr_thresh1)
      colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "meanLogFC_FDR_co2"] <- paste0("meanLogFC_thresh_FDR",
                                                                                            fdr_thresh2)
      
      
      colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "meanCorr_FDR_co1"] <- paste0("meanCorr_thresh_FDR",
                                                                                            fdr_thresh1)
      colnames(exprds_hicds_dt)[colnames(exprds_hicds_dt) == "meanCorr_FDR_co2"] <- paste0("meanCorr_thresh_FDR",
                                                                                            fdr_thresh2)
      
      
      exprds_hicds_dt <- exprds_hicds_dt[order(exprds_hicds_dt$adjPvalComb),]
      
      exprds_hicds_dt
    } # end-foreach iterating over exprds
    hicds_dt
  } # end-foreach iterating over hicds
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  save(all_result_dt, file = outFile)
  cat(paste0("... written: ", outFile,  "\n"))
  
  
  takecols <- colnames(all_result_dt)[! colnames(all_result_dt) %in% c(paste0("meanLogFC_thresh_FDR",
                                                                              fdr_thresh1),
                                                                       paste0("meanLogFC_thresh_FDR",
                                                                              fdr_thresh2),
                                                                       paste0("meanCorr_thresh_FDR",
                                                                              fdr_thresh1),
                                                                       paste0("meanCorr_thresh_FDR",
                                                                              fdr_thresh2))]
  
  
  all_result_dt_txt <- all_result_dt[, takecols]
  outFile <- file.path(outFolder, "all_result_dt.txt")
  write.table(all_result_dt_txt, file = outFile, sep="\t", quote=F, row.names=F, col.names=T,append=F)
  cat(paste0("... written: ", outFile,  "\n"))

    
} else {
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  all_result_dt <- eval(parse(text = load(outFile)))
}

############################################################################################# PLOT N SIGNIF

pvalThresh005 <- 0.05
pvalThresh001 <- 0.01

idvars <- c("hicds", "exprds", "dataset")

all_result_dt$dataset <- paste0(all_result_dt$hicds, "\n", all_result_dt$exprds)

all_result_dt$signifFDR02_adjPvalComb001 <- all_result_dt$adjPvalComb<=pvalThresh001 & all_result_dt$signifFDR_0.2

nSignif_empPval005 <- aggregate(adjPvalComb ~ dataset+hicds+exprds, FUN=function(x) sum(x<=pvalThresh005), data=all_result_dt)
colnames(nSignif_empPval005)[colnames(nSignif_empPval005) == "adjPvalComb"] <- "adjPvalComb005"

nSignif_empPval001 <- aggregate(adjPvalComb ~ dataset+hicds+exprds, FUN=function(x) sum(x<=pvalThresh001), data=all_result_dt)
colnames(nSignif_empPval001)[colnames(nSignif_empPval001) == "adjPvalComb"] <- "adjPvalComb001"

nSignif_fdr02_empPval001 <-  aggregate(signifFDR02_adjPvalComb001 ~ dataset+hicds+exprds, data=all_result_dt,  sum)
nSignif_fdr01 <-  aggregate(signifFDR_0.1 ~ dataset+hicds+exprds, data=all_result_dt,  sum)
nSignif_fdr02 <-  aggregate(signifFDR_0.2 ~ dataset+hicds+exprds, data=all_result_dt,  sum)

nSignif_dt <- merge(nSignif_fdr02_empPval001, 
                    merge(nSignif_empPval001, 
                          merge(nSignif_empPval005, merge(nSignif_fdr02, nSignif_fdr01, by=idvars, all=TRUE),
                                by =idvars, all=TRUE), 
                          by =idvars, all=TRUE),
                    by =idvars, all=TRUE)
nSignif_dt <- nSignif_dt[order(nSignif_dt$adjPvalComb005, decreasing = TRUE),]
ds_levels <- as.character(nSignif_dt$dataset)

all_dt <- melt(nSignif_dt, id=idvars)

all_dt$dataset <- factor(all_dt$dataset, levels=ds_levels)
all_dt <- all_dt[order(as.numeric(all_dt$dataset)),]
all_dt$exprds_type_col <- all_cols[all_cmps[all_dt$exprds]]
mycols <- all_dt$exprds_type_col[as.character(all_dt$variable) == "adjPvalComb005"] 

all_dt$variable <- factor(all_dt$variable, levels=c("adjPvalComb001", "adjPvalComb005", "signifFDR_0.1", "signifFDR_0.2", "signifFDR02_adjPvalComb001"))
stopifnot(!is.na(all_dt$variable))
############################################################
### => all variables; adjPvalComb001 sorted
############################################################
p_var <-  ggplot(all_dt, aes(x = dataset, y = value, fill = variable)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE) +
  ggtitle("# signif. TADs", subtitle = "(adj. emp. p-val. sorted)")+
  scale_x_discrete(name="")+
  labs(fill="")+
  scale_fill_manual(values=c(adjPvalComb001=pval001_col, adjPvalComb005=pval005_col, signifFDR_0.1=fdr01_col, signifFDR_0.2=fdr02_col, signifFDR02_adjPvalComb001=fdr_pval_col), 
                    labels=c("signif. adjPvalComb<=0.001","signif. adjPvalComb<=0.005", "signif. FDR<=0.1", "signif. FDR<=0.2", "signif. FDR<=0.2, adjPvalComb<=0.01"))+
  scale_y_continuous(name=paste0("# signif. TADs"),
                     breaks = scales::pretty_breaks(n = 10))+
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

outFile <- file.path(outFolder, paste0("all_ds_nSignif_empPval_sorted.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

############################################################
### => only FDR; FDR 0.2 sorted
############################################################
nSignif_dt <- nSignif_dt[order(nSignif_dt$signifFDR_0.2, decreasing = TRUE),]
ds_levels_fdr <- as.character(nSignif_dt$dataset)

all_dt$dataset <- factor(all_dt$dataset, levels=ds_levels_fdr)
all_dt <- all_dt[order(as.numeric(all_dt$dataset)),]
all_dt$exprds_type_col <- all_cols[all_cmps[all_dt$exprds]]
mycols <- all_dt$exprds_type_col[as.character(all_dt$variable) == "signifFDR_0.2"]

p_var <-  ggplot(all_dt[!grepl( "adjPvalComb", all_dt$variable),], aes(x = dataset, y = value, fill = variable)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE) +
  ggtitle("# signif. TADs", subtitle = "(FDR 0.2 sorted)")+
  scale_x_discrete(name="")+
  labs(fill="")+
  scale_fill_manual(values=c(signifFDR_0.1=fdr01_col, signifFDR_0.2=fdr02_col), 
                    labels=c("signif. FDR<=0.1", "signif. FDR<=0.2"))+
  scale_y_continuous(name=paste0("# signif. TADs"),
                     breaks = scales::pretty_breaks(n = 10))+
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

outFile <- file.path(outFolder, paste0("all_ds_nSignif_fdr0.2_sorted.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))
############################################################
### => onyl FDR;  FDR 0.1 sorted
############################################################

nSignif_dt <- nSignif_dt[order(nSignif_dt$signifFDR_0.1, decreasing = TRUE),]
ds_levels_fdr <- as.character(nSignif_dt$dataset)

all_dt$dataset <- factor(all_dt$dataset, levels=ds_levels_fdr)
all_dt <- all_dt[order(as.numeric(all_dt$dataset)),]
all_dt$exprds_type_col <- all_cols[all_cmps[all_dt$exprds]]
mycols <- all_dt$exprds_type_col[as.character(all_dt$variable) == "signifFDR_0.1"]


p_var <-  ggplot(all_dt[!grepl( "adjPvalComb", all_dt$variable),], aes(x = dataset, y = value, fill = variable)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE) +
  ggtitle("# signif. TADs", subtitle = "(FDR 0.1 sorted)")+
  scale_x_discrete(name="")+
  labs(fill="")+
  scale_fill_manual(values=c(signifFDR_0.1=fdr01_col, signifFDR_0.2=fdr02_col), 
                    labels=c("signif. FDR<=0.1", "signif. FDR<=0.2"))+
  scale_y_continuous(name=paste0("# signif. TADs"),
                     breaks = scales::pretty_breaks(n = 10))+
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

outFile <- file.path(outFolder, paste0("all_ds_nSignif_fdr0.1_sorted.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

############################################################
### => onyl EMPPVAL;   0.05 sorted
############################################################

nSignif_dt <- nSignif_dt[order(nSignif_dt$adjPvalComb005, decreasing = TRUE),]
ds_levels_fdr <- as.character(nSignif_dt$dataset)

all_dt$dataset <- factor(all_dt$dataset, levels=ds_levels_fdr)
all_dt <- all_dt[order(as.numeric(all_dt$dataset)),]
all_dt$exprds_type_col <- all_cols[all_cmps[all_dt$exprds]]
mycols <- all_dt$exprds_type_col[as.character(all_dt$variable) == "adjPvalComb005"]


p_var <-  ggplot(all_dt[!grepl( "FDR", all_dt$variable),], aes(x = dataset, y = value, fill = variable)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE) +
  ggtitle("# signif. TADs", subtitle = "(empPval 0.05 sorted)")+
  scale_x_discrete(name="")+
  labs(fill="")+
  scale_fill_manual(values=c(adjPvalComb001=pval001_col, adjPvalComb005=pval005_col), 
                    labels=c("signif. adjPvalComb<=0.001","signif. adjPvalComb<=0.005"))+
  scale_y_continuous(name=paste0("# signif. TADs"),
                     breaks = scales::pretty_breaks(n = 10))+
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

outFile <- file.path(outFolder, paste0("all_ds_nSignif_empPval0.01_empPval0.05_sorted.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


############################################################
### => only pval0.001 + FDR 0.2, sorted
############################################################



nSignif_dt <- nSignif_dt[order(nSignif_dt$signifFDR02_adjPvalComb001, decreasing = TRUE),]
ds_levels_fdr <- as.character(nSignif_dt$dataset)

all_dt$dataset <- factor(all_dt$dataset, levels=ds_levels_fdr)
all_dt <- all_dt[order(as.numeric(all_dt$dataset)),]
all_dt$exprds_type_col <- all_cols[all_cmps[all_dt$exprds]]
mycols <- all_dt$exprds_type_col[as.character(all_dt$variable) == "signifFDR02_adjPvalComb001"]

p_var <-  ggplot(all_dt[all_dt$variable == "signifFDR02_adjPvalComb001",], aes(x = dataset, y = value, fill = variable)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE) +
  ggtitle("# signif. TADs", subtitle = "(FDR 0.2 + emp. p-val 0.001 sorted)")+
  scale_x_discrete(name="")+
  labs(fill="")+
  scale_fill_manual(values=c("signifFDR02_adjPvalComb001"=fdr_pval_col), 
                    labels=c("signif. FDR<=0.2, adjPvalComb<=0.01"))+
  scale_y_continuous(name=paste0("# signif. TADs"),
                     breaks = scales::pretty_breaks(n = 10))+
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

outFile <- file.path(outFolder, paste0("all_ds_nSignif_empPval0.001_fdr0.2_sorted.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))





#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))










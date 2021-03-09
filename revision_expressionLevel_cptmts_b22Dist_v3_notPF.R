require(ggsci)
require(ggpubr)
require(ggplot2)
require(ggrepel)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("revision_settings.R")
require(stringr)
require(foreach)
require(doMC)
registerDoMC(40)
plotType <- "png"
myWidthGG <- 7
myHeightGG <- 5
myHeight <- 400
myWidth <- 400

set.seed(020321)

plotCex <- 1.2

tadSignifThresh <- 0.01

### SAMPLE SAME SIZE DIST AS IN B22

nRandom <- 1000
runPermut <- TRUE
# Rscript revision_expressionLevel_cptmts_b22Dist_v3_notPF.R

outFolder <- file.path("REVISION_EXPRESSIONLEVEL_CPTMTS_B22DIST_V3_NOTPF")
dir.create(outFolder, recursive=TRUE)

#### #### #### #### #### #### #### #### #### #### #### #### 
#### retrieve the purity tagged tads

runFolder <- "."

purity_ds <- "aran"
pm <- "CPE"
purity_plot_name <- paste0("Aran - ", pm)

corMet <- "pearson"
transfExpr <- "log10"
corrPurityQtThresh <- 0.05
signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)


purity_file <- file.path(runFolder,"ALLTADS_AND_PURITY_FINAL", purity_ds, pm, transfExpr, "all_ds_corrPurity_dt.Rdata") # here _final INPUT
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)
agg_purity$region_ID <- file.path(agg_purity$dataset, agg_purity$region)
stopifnot(!duplicated(agg_purity$region_ID))

result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!
merge_dt$region_ID <- file.path(merge_dt$dataset, merge_dt$region)
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", signifThresh), paste0("adj. p-val >", signifThresh) )
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))

# purity_flagged_tads <- merge_dt$region_ID[merge_dt$purityCorr <= purityCorrThresh]
tokeep_tads <- merge_dt$region_ID[merge_dt$purityCorr > purityCorrThresh]
cat(paste0( "purityCorrThresh = ", round(purityCorrThresh, 4), "\n"))
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 




inFolder <- "REVISION_CHANGES_CPTMTLABELS_ALLDS"
outFile <- file.path(inFolder, "tad2cptmt_dt.Rdata")
tad2cptmt_dt <- get(load(outFile))

stopifnot(tad2cptmt_dt$tad_binaryCptmtLab == tad2cptmt_dt$start_binaryCptmtLab)
stopifnot(tad2cptmt_dt$tad_eightCptmtLab == tad2cptmt_dt$start_eightCptmtLab)

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$region_ID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
stopifnot(!duplicated(final_table_DT$regionID))
regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)
signif_tads <- final_table_DT$regionID[final_table_DT$adjPvalComb <= tadSignifThresh]
final_table_DT$nGenes <- str_count(final_table_DT$region_genes, pattern=",") + 1
stopifnot(!is.na(final_table_DT$nGenes))
stopifnot(final_table_DT$nGenes >= 3)
final_table_DT$size_log10 <- log10(final_table_DT$end-final_table_DT$start+1)
stopifnot(!is.na(final_table_DT$size_log10))

ngenes_values <- setNames(final_table_DT$nGenes, final_table_DT$region_ID)
size_values <- setNames(final_table_DT$size_log10, final_table_DT$region_ID)

logFC_values <- setNames(final_table_DT$meanLogFC, final_table_DT$region_ID)
corr_values <- setNames(final_table_DT$meanCorr, final_table_DT$region_ID)

stopifnot(setequal(names(logFC_values), tad2cptmt_dt$region_ID))
stopifnot(setequal(names(corr_values), tad2cptmt_dt$region_ID))

tad2cptmt_dt[,paste0("meanCorr")] <- corr_values[tad2cptmt_dt$region_ID]
tad2cptmt_dt[,paste0("meanLogFC")] <- logFC_values[tad2cptmt_dt$region_ID]

tad2cptmt_dt[,paste0("sizeNgenes")] <- ngenes_values[tad2cptmt_dt$region_ID]
tad2cptmt_dt[,paste0("sizeBp")] <- size_values[tad2cptmt_dt$region_ID]


legTitle <- ""

tad2cptmt_dt$tad_binaryCptmtLab <- factor(tad2cptmt_dt$tad_binaryCptmtLab,
                                          levels=as.character(sort(unique(as.character(
                                            tad2cptmt_dt$tad_binaryCptmtLab)))))
stopifnot(!is.na(tad2cptmt_dt$tad_binaryCptmtLab))


tad2cptmt_dt$tad_eightCptmtLab <- factor(tad2cptmt_dt$tad_eightCptmtLab,
                                         levels=as.character(sort(unique(as.character(
                                           tad2cptmt_dt$tad_eightCptmtLab)))))
stopifnot(!is.na(tad2cptmt_dt$tad_eightCptmtLab))


########################################## ADD HERE THE PURTIY FILTER

stopifnot(tokeep_tads %in% tad2cptmt_dt$region_ID)
tad2cptmt_dt <- tad2cptmt_dt[tad2cptmt_dt$region_ID %in% tokeep_tads,]


###################
### PREPARE THE purity corr DATA
###################

purityCorr_values <- setNames(agg_purity$purityCorr, agg_purity$region_ID)
stopifnot(tad2cptmt_dt$region_ID %in% names(purityCorr_values))
tad2cptmt_dt[,paste0("meanPurityCorr")] <- purityCorr_values[tad2cptmt_dt$region_ID]


###################
### PREPARE THE GENE FC DATA
###################

expr_var <- "aggLog10Expr"
all_exprLevel_dt <- get(load(file.path("REVISION_EXPRESSION_LEVEL", paste0(expr_var, "_aggByTAD_mean.Rdata"))))
exprVar_values <- setNames(all_exprLevel_dt[,paste0(expr_var)], all_exprLevel_dt$regionID)
# stopifnot(setequal(names(exprVar_values), tad2cptmt_dt$region_ID)) # not true PF
stopifnot(tad2cptmt_dt$region_ID %in% names(exprVar_values)) # not true PF
tad2cptmt_dt[,paste0(expr_var)] <- exprVar_values[tad2cptmt_dt$region_ID]


############################################################################
############################################################################

# perform here the random sampling
stopifnot(!duplicated(tad2cptmt_dt$region_ID))

cptmt_var <- "tad_eightCptmtLab"
ref_cptmt <- "B.2.2"
strat_var <- "sizeNgenes"

tad2cptmt_dt[,cptmt_var] <- as.character(tad2cptmt_dt[,cptmt_var])
stopifnot(ref_cptmt %in% tad2cptmt_dt[,cptmt_var] )

ref_sizeDist <-table( tad2cptmt_dt[,strat_var][tad2cptmt_dt[,cptmt_var] == ref_cptmt])

print("ref_sizeDist\n")
print(ref_sizeDist)
print("sum ref_sizeDist\n")
print(sum(ref_sizeDist))
stop("-ok \n")

stopifnot(sum(ref_sizeDist) == sum(tad2cptmt_dt[,cptmt_var] == ref_cptmt))

othercptmt_dt <- tad2cptmt_dt[tad2cptmt_dt[,cptmt_var] != ref_cptmt,]

nOtherCptmts <- length(unique(othercptmt_dt[,cptmt_var]))

refcptmt_dt <- tad2cptmt_dt[tad2cptmt_dt[,cptmt_var] == ref_cptmt,]

stopifnot(!ref_cptmt %in% othercptmt_dt[,cptmt_var])

cat(paste0("... start permut\n"))

if(runPermut) {
  rd_tad2cptmt_dt <- foreach(i = 1:nRandom, .combine='rbind') %dopar% {
    # for each other compartment sample according to size dist
    sub_dt=othercptmt_dt[othercptmt_dt[,cptmt_var] == othercptmt_dt[1,cptmt_var] ,]
    rd_dt <- do.call(rbind, by(othercptmt_dt,othercptmt_dt[,cptmt_var], function(sub_dt) {
      
      stopifnot(length(unique(sub_dt[,cptmt_var] )) == 1)
      
      i_size=1
      cptmt_dt <- do.call(rbind, lapply(seq_along(ref_sizeDist), function(i_size) {
        
        tadsize <- as.numeric(as.character(names(ref_sizeDist)[i_size]))
        tosample <- as.numeric(ref_sizeDist[i_size])
        size_sub_dt <- sub_dt[sub_dt[,strat_var] == tadsize,]    
        rownames(size_sub_dt) <- NULL
        stopifnot(nrow(size_sub_dt) > 0)
        
        sampleReplace <- (tosample > nrow(size_sub_dt)) # sample with replacement only if more to sample actual # of rows
        
        sampled_rows <- sample(1:nrow(size_sub_dt), size = tosample, replace=sampleReplace)
        out_dt <- size_sub_dt[sampled_rows,]
        stopifnot(nrow(out_dt) == tosample)
        stopifnot(out_dt[,strat_var] == tadsize)
        out_dt
      })) 
      stopifnot(sum(cptmt_dt[,strat_var]) == sum(refcptmt_dt[,strat_var]))
      stopifnot(table(cptmt_dt[,strat_var]) == table(refcptmt_dt[,strat_var]))
      cptmt_dt
    } ))
    stopifnot(nrow(rd_dt) == nOtherCptmts * sum(ref_sizeDist))
    rd_dt
  }
  stopifnot(nrow(rd_tad2cptmt_dt) == nOtherCptmts * nRandom * sum(ref_sizeDist))
  outFile <- file.path(outFolder, "rd_tad2cptmt_dt.Rdata")
  # save(rd_tad2cptmt_dt, file = outFile, version=2)
} else {
  outFile <- file.path(outFolder, "rd_tad2cptmt_dt.Rdata")
  rd_tad2cptmt_dt <- get(load(outFile))
}

cat(paste0("... end permut\n"))

pseudo_tad2cptmt_dt <- rbind(rd_tad2cptmt_dt, refcptmt_dt)

tad2cptmt_dt <- pseudo_tad2cptmt_dt

# for each of the other compartments, sample same amount of dist size


plotTit <- paste0(cptmt_var, " stratif. by ", strat_var, " according to ", ref_cptmt)

all_cptmt_vars <- c("tad_binaryCptmtLab","tad_eightCptmtLab")

all_plot_vars <- c(expr_var, "meanCorr", "meanLogFC", "sizeNgenes", "sizeBp",  "meanPurityCorr")


for(cptmt_var in all_cptmt_vars){
  
  for(plot_var in all_plot_vars) {
    
    # plotTit <- paste0(plot_var, " by ", cptmt_var)
    
    mySub <- paste0("# DS = ", length(unique(file.path(tad2cptmt_dt$hicds,tad2cptmt_dt$exprds))),
                    "; # TADs = ", nrow(tad2cptmt_dt), "; # permut. = ", nRandom)
    
    pbox <- ggboxplot(tad2cptmt_dt, 
                      outlier.shape=NA,
                      x=paste0(cptmt_var),
                      y=paste0(plot_var),
                      add="jitter") + 
      mytheme +
      ggtitle(plotTit, subtitle = mySub)+
      # scale_color_manual(values=my_cols)+
      # scale_fill_manual(values=my_cols)  +
      labs(color=paste0(legTitle),fill=paste0(legTitle), x="", y=paste0("TAD ", plot_var))+
      # guides(color=FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
    
    outFile <- file.path(outFolder,paste0(plot_var, "_byCptmt_", cptmt_var, "_boxplot.", plotType))
    ggsave(pbox, filename = outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    mySub <- paste0("# DS = ", length(unique(file.path(tad2cptmt_dt$hicds,tad2cptmt_dt$exprds))),
                    "; # TADs = ", nrow(tad2cptmt_dt), "; # permut. = ", nRandom)
    
    pbox <- ggviolin(tad2cptmt_dt, 
                     outlier.shape=NA,
                     add = "boxplot",
                     x=paste0(cptmt_var),
                     y=paste0(plot_var) #,
                     # add="jitter"
    ) + 
      mytheme +
      ggtitle(plotTit, subtitle = mySub)+
      # scale_color_manual(values=my_cols)+
      # scale_fill_manual(values=my_cols)  +
      labs(color=paste0(legTitle),fill=paste0(legTitle), x="", y=paste0("TAD ", plot_var))+
      # guides(color=FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
    
    outFile <- file.path(outFolder,paste0(plot_var, "_byCptmt_", cptmt_var, "_violinplot.", plotType))
    ggsave(pbox, filename = outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    if(cptmt_var == "tad_eightCptmtLab" ){
      
      tad2cptmt_dt$cptmt_var_grouped <- gsub("(.).+", "\\1", tad2cptmt_dt$tad_eightCptmtLab)
      tad2cptmt_dt$cptmt_var_grouped[as.character(tad2cptmt_dt$tad_eightCptmtLab) == "B.2.2"] <- "B.2.2"
      
      # save(tad2cptmt_dt, file="tmp_tad2cptmt_dt.Rdata", version=2)
      
      p3 <- ggdensity(tad2cptmt_dt,
                      x = paste0(plot_var),
                      # y = "..count..",
                      y = "..density..",
                      ylab="Density",
                      rug = FALSE,                      # Add marginal rug
                      color = paste0("cptmt_var_grouped"),
                      fill = paste0("cptmt_var_grouped"),
                      palette = "jco"
      ) +
        labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density", x=paste0("TAD ", plot_var))+
        
        ggtitle(plotTit, subtitle = mySub)+
        guides(color=FALSE)+
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
        scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
        mytheme
      
      outFile <- file.path(outFolder,paste0(plot_var, "_byCptmt_", cptmt_var, "_groupedDensityplot.", plotType))
      ggsave(p3, filename = outFile, height=myHeightGG, width=myWidthGG)
      cat(paste0("... written: ", outFile, "\n"))
    }
    
    
    
  }
}

# ############## DENSPLOT WITH THE CONTINUOUS RANKS
# 
# all_cptmt_vars <- c("tadCptmtNormRank")
# all_plot_vars <- c(expr_var, "meanCorr", "meanLogFC", "sizeNgenes", "sizeBp")
# 
# cptmt_var="tadCptmtNormRank"
# plot_var="meanCorr"
# 
# for(cptmt_var in all_cptmt_vars){
#   
#   for(plot_var in all_plot_vars) {
#     
#     
#     myx <- tad2cptmt_dt[,paste0(cptmt_var)]
#     myy <- tad2cptmt_dt[,paste0(plot_var)]
# 
#         
#     # plotTit <- paste0( plot_var, " by ", cptmt_var)
#     
#     mySub <- paste0("# DS = ", length(unique(file.path(tad2cptmt_dt$hicds,tad2cptmt_dt$exprds))),
#                     "; # TADs = ", nrow(tad2cptmt_dt))
#     
#     outFile <- file.path(outFolder,paste0(plot_var, "_byCptmt_", cptmt_var, "_densplot.", plotType))
#     do.call(plotType, list(outFile, height=myWidth, width=myWidth))
#     densplot(x=myx,
#              y=myy,
#              xlab = paste0(cptmt_var),
#              ylab = paste0(plot_var),
#              main  = plotTit,
#              plot.main=plotCex,
#              plot.main=plotCex,
#              plot.main=plotCex
#     )
#     mtext(side=3, text=mySub)
#     addCorr(x=myx,y=myy, legPos="topleft", bty="n", corMet="spearman")
#     foo <- dev.off()
#     cat(paste0("... written: ", outFile, "\n"))
#     
# 
#   }
# }
# 

# ##### pie chart signif by cptmts
# 
# stopifnot(!duplicated(tad2cptmt_dt$region_ID))
# 
# signif_dt <- tad2cptmt_dt[tad2cptmt_dt$adjPvalComb <= tadSignifThresh,]
# 
# all_cptmt_vars <- c("tad_binaryCptmtLab","tad_eightCptmtLab")
# 
# cptmt_var="tad_binaryCptmtLab"
# 
# 
# for(cptmt_var in all_cptmt_vars){
#   
#   tmp_dt <- aggregate(as.formula(paste0("region_ID ~ ", cptmt_var)), data=tad2cptmt_dt, FUN=length)
#   colnames(tmp_dt)[colnames(tmp_dt) == "region_ID"] <- "nTADs"
#   tmp_dt$ratioTADs <- tmp_dt$nTADs/nrow(tad2cptmt_dt)
#   tmp_dt$ratioTADs_lab <- paste0(round(tmp_dt$ratioTADs*100, 2), "%")
#   
#   myTit <- paste0("Dist. of all TADs across ", cptmt_var)
#   mysub <- paste0("# DS = ", length(unique(file.path(tad2cptmt_dt$hicds, tad2cptmt_dt$exprds))) , "; # TADs = ", nrow(tad2cptmt_dt))
#   
#   p_dist <- ggplot(tmp_dt, aes_string(x="1", y="ratioTADs", fill=cptmt_var)) +
#     geom_col() +
#     geom_text_repel(aes(label = ratioTADs_lab), position = position_stack(vjust = 0.5))+
#     coord_polar(theta = "y") + 
#     ggtitle(myTit, subtitle=mysub) +
#     theme_void() +    
#     labs(fill="") +
#     blank_theme
#   
#   outFile <- file.path(outFolder, paste0( "distAllTADs_by_", cptmt_var, "_pie.", plotType))
#   ggsave(p_dist, file=outFile, height=myHeightGG, width=myWidthGG)
#   cat(paste0("... written: ", outFile, "\n"))
#   
#   
#   agg_dt <- aggregate(as.formula(paste0("region_ID ~ ", cptmt_var)), data=signif_dt, FUN=length)
#   colnames(agg_dt)[colnames(agg_dt) == "region_ID"] <- "nTADs"
#   agg_dt$ratioTADs <- agg_dt$nTADs/nrow(signif_dt)
#   
#   agg_dt$ratioTADs_lab <- paste0(round(agg_dt$ratioTADs*100, 2), "%")
#   
#   myTit <- paste0("Dist. signif. TADs across ", cptmt_var)
#   mysub <- paste0("# DS = ", length(unique(file.path(signif_dt$hicds, signif_dt$exprds))) , "; # TADs = ", nrow(signif_dt))
#   
#   p_dist <- ggplot(agg_dt, aes_string(x="1", y="ratioTADs", fill=cptmt_var)) +
#     geom_col() +
#     geom_text_repel(aes(label = ratioTADs_lab), position = position_stack(vjust = 0.5))+
#     coord_polar(theta = "y") + 
#     ggtitle(myTit, subtitle=mysub) +
#     theme_void() +    
#     labs(fill="") +
#     blank_theme
#   
#   outFile <- file.path(outFolder, paste0( "distSignifTADs_by_", cptmt_var, "_pie.", plotType))
#   ggsave(p_dist, file=outFile, height=myHeightGG, width=myWidthGG)
#   cat(paste0("... written: ", outFile, "\n"))
#   
# }
# 
# 
# 
# 
# # hicds_norm ="LI_40kb"
# # 
# # all_norm_files <- list.files(file.path(hicds_norm, "FINAL_DOMAINS_WITH_SCORES"), pattern="final_domains_with_scores.txt$", full.names = TRUE)
# normFile = all_norm_files[1]
# dt <- read.delim(normFile, header=F, col.names=c("chromo", "start", "end", "rankValue"))
# 
# x=tad2cptmt_final_dt[tad2cptmt_final_dt$hicds == "LI_40kb" & grepl("chr1_TAD", tad2cptmt_final_dt$region),]
# 
# yy = merge(x, dt, by=c("start", "end"))
# 
# 
# # aggregate the rank values
# allChr_norm_rankDT <- foreach(normFile = all_norm_files, .combine='rbind') %dopar% {
#   dt <- read.delim(normFile, header=F, col.names=c("chromo", "start", "end", "rankValue"))
#   dt
# }
# > sum(abs(yy$rankValue - yy$startCptmtNormRank)
#       + )
# [1] 0




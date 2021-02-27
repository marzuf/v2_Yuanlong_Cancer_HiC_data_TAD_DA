# Rscript revision_expressionLevel_cptmts.R

inFolder <- "REVISION_CHANGES_CPTMTLABELS_ALLDS"
outFile <- file.path(inFolder, "tad2cptmt_dt.Rdata")
tad2cptmt_dt <- get(load(outFile))

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$region_ID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
stopifnot(!duplicated(final_table_DT$regionID))
regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)
signif_tads <- final_table_DT$regionID[final_table_DT$adjPvalComb <= tadSignifThresh]

logFC_values <- setNames(final_table_DT$meanLogFC, final_table_DT$region_ID)
corr_values <- setNames(final_table_DT$meanCorr, final_table_DT$region_ID)

stopifnot(setequal(names(logFC_values), tad2cptmt_dt$region_ID))
stopifnot(setequal(names(corr_values), tad2cptmt_dt$region_ID))

tad2cptmt_dt[,paste0("corr_values")] <- corr_values[tad2cptmt_dt$region_ID]
tad2cptmt_dt[,paste0("logFC_values")] <- logFC_values[tad2cptmt_dt$region_ID]

ggboxplot(tad2cptmt_dt, 
          x="tad_binaryCptmtLab",
          y=paste0("corr_values")) + 
  mytheme

ggboxplot(tad2cptmt_dt, 
          x="tad_binaryCptmtLab",
          y=paste0("logFC_values")) + 
  mytheme

###################
### PREPARE THE GENE FC DATA
###################

expr_var <- "aggLog10Expr"
all_exprLevel_dt <- get(load(file.path("REVISION_EXPRESSION_LEVEL", paste0(expr_var, "_aggByTAD_mean.Rdata"))))
exprVar_values <- setNames(all_exprLevel_dt[,paste0(expr_var)], all_exprLevel_dt$regionID)
stopifnot(setequal(names(exprVar_values), tad2cptmt_dt$region_ID))
tad2cptmt_dt[,paste0(expr_var)] <- exprVar_values[tad2cptmt_dt$region_ID]


plotTit <- ""
mysub <- ""

# set correct levels
tad2cptmt_dt$tad_eightCptmtLab <- factor(tad2cptmt_dt$tad_eightCptmtLab,
                                         levels=sort(as.character(unique(tad2cptmt_dt$tad_eightCptmtLab))))
stopifnot(!is.na(tad2cptmt_dt$tad_eightCptmtLab))
tad2cptmt_dt$tad_binaryCptmtLab <- factor(tad2cptmt_dt$tad_binaryCptmtLab,
                                         levels=sort(as.character(unique(tad2cptmt_dt$tad_binaryCptmtLab))))
stopifnot(!is.na(tad2cptmt_dt$tad_binaryCptmtLab))

all_cptmt_vars <- c("tad_binaryCptmtLab","tad_eightCptmtLab")

all_plot_vars <- c(expr_var, "corr_values", "logFC_values")

cptmt_var="tad_eightCptmtLab"
plot_var="corr_values"
plot_var=expr_var
for(cptmt_var in all_cptmt_vars){
  
  for(plot_var in all_plot_vars) {
    
    
    ggboxplot(tad2cptmt_dt, 
              x=paste0(cptmt_var),
              y=paste0(plot_var)) + 
      mytheme
    
    
    outFile <- file.path(outFolder, paste0("tad_genes_annot_", plot_var, "_signif_notsignif_density.", plotType))
    ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
  }
}




all_cptmt_vars <- c("tadCptmtNormRank")

all_plot_vars <- c(expr_var, "corr_values", "logFC_values")

cptmt_var="tadCptmtNormRank"
plot_var="corr_values"

for(cptmt_var in all_cptmt_vars){
  
  for(plot_var in all_plot_vars) {
    
    
    myx <- tad2cptmt_dt[,paste0(cptmt_var)]
    myy <- tad2cptmt_dt[,paste0(plot_var)]
    
    densplot(
              x=myx,
              y=myy)
                
    
    

    
    
  }
}





# hicds_norm ="LI_40kb"
# 
# all_norm_files <- list.files(file.path(hicds_norm, "FINAL_DOMAINS_WITH_SCORES"), pattern="final_domains_with_scores.txt$", full.names = TRUE)
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




# FIGURE 2: scatterplot meanFC meanCorr color by pval

options(scipen=100)

setDir = ""

# Rscript report_figure2.R   # 

script_name <- "report_figure2.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(ggplot2)
require(ggpubr)
require(ggforce)
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

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- ifelse(plotType=="png", 600, 10)
plotCex <- 1.2

myWidthGG <- 12
myHeightGG <- 12

outFolder <- "REPORT_FIGURE2"
dir.create(outFolder, recursive=TRUE)

signifTAD_thresh <- 0.05
signifTAD_FDR <- 0.2

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

final_table_DT$pval_signif <- final_table_DT$adjPvalComb <= signifTAD_thresh
final_table_DT$fdr_signif <- final_table_DT[, paste0("signifFDR_", signifTAD_FDR)]
final_table_DT$pval_fdr_signif <- final_table_DT$fdr_signif & final_table_DT$pval_signif

final_table_DT$adjPvalComb_log10 <- -log10(final_table_DT$adjPvalComb)

format_my_plot <- function(p) {
  p <- p + 
    gradient_color("red") + 
    # gradient_color("YlOrRd") + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    theme( 
      plot.subtitle = element_text(hjust=0.5),
      axis.text.y = element_text(color="black", size=12),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.text.x = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=14, face="bold")
    )
  return(p)
}

format_my_plot_facet <- function(p) {
  p <- format_my_plot(p)
  p <- p + theme(strip.text.x = element_text(size = 14))
  return(p)
}



p1 <- ggscatter(final_table_DT, 
                subtitle=paste0("n=", nrow(final_table_DT)),
               x = "meanLogFC", 
               y = "meanCorr",
               color = "adjPvalComb_log10")
p1 <- format_my_plot(p1)
  
p1b <- ggscatter(final_table_DT[final_table_DT$adjPvalComb <=0.01,], 
                 subtitle=paste0("n=", nrow(final_table_DT[final_table_DT$adjPvalComb <=0.01,])),
               x = "meanLogFC", 
               y = "meanCorr",
               color = "adjPvalComb_log10")
p1b <- format_my_plot(p1b)

final_table_DT$cmpType <- all_cmps[final_table_DT$exprds]
stopifnot(!is.na(final_table_DT$cmpType))

p2 <- ggscatter(final_table_DT, 
          x = "meanLogFC",
          y = "meanCorr", 
          color = "adjPvalComb_log10",
          scales="free",
          facet.by = "cmpType")
p2 <- format_my_plot_facet(p2)


for(cmp_type in unique(as.character(final_table_DT$cmpType))) {
  p2annot_1 <-  data.frame(meanLogFC = min(final_table_DT$meanLogFC[final_table_DT$cmpType == cmp_type]),
                           meanCorr = min(final_table_DT$meanCorr[final_table_DT$cmpType == cmp_type]),
                           lab = paste0("n=", sum(final_table_DT$cmpType == cmp_type), 
                           "; #DS=", length(unique(paste0(final_table_DT[final_table_DT$cmpType == cmp_type,]$hicds, final_table_DT[final_table_DT$cmpType == cmp_type,]$exprds)))),
                           cmpType = paste0(cmp_type))
  p2 <- p2 + geom_text(data = p2annot_1, label = p2annot_1$lab, hjust=0, vjust=0)
}



tmpdt <- final_table_DT[final_table_DT$adjPvalComb <=0.01,]

p2b <- ggscatter(tmpdt, 
                x = "meanLogFC",
                y = "meanCorr", 
                color = "adjPvalComb_log10",
                scales="free",
                facet.by = "cmpType")
p2b <- format_my_plot_facet(p2b)

for(cmp_type in unique(as.character(tmpdt$cmpType))) {
  p2bannot_1 <-  data.frame(meanLogFC = min(tmpdt$meanLogFC[tmpdt$cmpType == cmp_type]),
                           meanCorr = min(tmpdt$meanCorr[tmpdt$cmpType == cmp_type]),
                           lab = paste0("n=", sum(tmpdt$cmpType == cmp_type), "; #DS=", length(unique(paste0(tmpdt[tmpdt$cmpType == cmp_type,]$hicds, tmpdt[tmpdt$cmpType == cmp_type,]$exprds)))),
                           cmpType = paste0(cmp_type))
  p2b <- p2b + geom_text(data = p2bannot_1, label = p2bannot_1$lab, hjust=0, vjust=0)
}



arr1 <- ggarrange(p1, p1b, ncol = 2, labels = c("ALL DATA", "ONLY SIGNIF. TADs")) #
arr2 <- ggarrange(p2, p2b, nrow = 2, labels = c("ALL DATA", "ONLY SIGNIF. TADs")) #


arr_all <- ggarrange(arr1, arr2,
          nrow = 2,
          heights=c(1/3,2/3)
          )

outFile <- file.path(outFolder, paste0("meanCorr_meanFC_allPlots.", plotType))
ggsave(plot = arr_all, filename = outFile, height=myHeightGG*2, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))


# density
signif_thresh <- 0.01

signif_dt <- final_table_DT[final_table_DT$adjPvalComb_log10 >= -log10(signif_thresh),]

subtit1 <- paste0(paste0("# ", names(table(final_table_DT$cmpType))), "=", as.numeric(table(final_table_DT$cmpType)), collapse="; ")

meanfc_p <- ggdensity(signif_dt,
          title=paste0("Mean logFC of signif. TADs"),
          subtitle=paste0("signif. thresh = ", signif_thresh, "\n", subtit1 ),
          x = "meanLogFC",
          # add = "mean", rug = TRUE,
          color = "cmpType", fill = "cmpType",
          palette = all_cols) + labs(fill="", color="")+
  theme(plot.title = element_text(face="bold", size=14))
outFile <- file.path(outFolder, paste0("signifTADs_meanFC_density.", plotType))
ggsave(plot = meanfc_p, filename = outFile, height=myHeightGG/2, width = myWidthGG/1.8)
cat(paste0("... written: ", outFile, "\n"))


meancorr_p <- ggdensity(signif_dt,
                      title=paste0("Mean corr. of signif. TADs"),
                      subtitle=paste0("signif. thresh = ", signif_thresh, "\n", subtit1 ),
                      x = "meanCorr",
                      # add = "mean", rug = TRUE,
                      color = "cmpType", fill = "cmpType",
                      palette = all_cols) + labs(fill="", color="")+
  theme(plot.title = element_text(face="bold", size=14))
outFile <- file.path(outFolder, paste0("signifTADs_meanCorr_density.", plotType))
ggsave(plot = meancorr_p, filename = outFile, height=myHeightGG/2, width = myWidthGG/1.8)
cat(paste0("... written: ", outFile, "\n"))

##########################################################################################
##########################################################################################
##########################################################################################

cat("*** DONE - ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))









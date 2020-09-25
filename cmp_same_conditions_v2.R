options(scipen=100)

setDir <- ""

# Rscript cmp_same_conditions_v2.R 

script_name <- "cmp_same_conditions_v2.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
require(reshape2)
require(ggpubr)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight

axisCex <- labCex <- 1.2

buildData <- TRUE

corMeth <- "pearson"

pipFolder <- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "CMP_SAME_CONDITIONS_V2"
dir.create(outFolder, recursive=TRUE)

gr_tr_dt <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))

all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[!( grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))

# retrieve the "duplicated" all_exprds
all_exprds_vect <- unlist(all_exprds)
dup_exprds <- unique(all_exprds_vect[duplicated(all_exprds_vect)])

expr_hicds <- sapply(dup_exprds, function(exprds) {
  names(all_exprds)[sapply(all_exprds, function(x) exprds %in% x )]
})

stopifnot(dup_exprds %in% names(expr_hicds))
  
### BUILD SIGNIF ALONG FDR THRESH
exprds = dup_exprds[1]
if(buildData){

  all_conds_dt <- foreach(exprds = dup_exprds, .combine='rbind') %dopar% {
    
    all_hicds <- expr_hicds[[paste0(exprds)]]
    
    all_hicds_comb <- combn(all_hicds,2)
    stopifnot(nrow(all_hicds_comb) == 2)
    stopifnot(ncol(all_hicds_comb) >= 1)
    i_cmb = 1
    exprds_dt <- foreach(i_cmb = 1:ncol(all_hicds_comb), .combine='rbind') %dopar% {
      
      
      hicds1 <- all_hicds_comb[1,i_cmb]
      hicds2 <- all_hicds_comb[2,i_cmb]
      
      gt1_dt <- gr_tr_dt[gr_tr_dt$hicds == hicds1 & gr_tr_dt$exprds == exprds,]
      
      # short check
      x_dt <- unique(gt1_dt[,c("region", "tad_adjCombPval", "tad_rank")])
      check_vect <- rank(x_dt$tad_adjCombPval, ties.method="min")
      stopifnot(check_vect == x_dt$tad_rank)
      
      
      gt2_dt <- gr_tr_dt[gr_tr_dt$hicds == hicds2 & gr_tr_dt$exprds == exprds,]
      
      merge_dt <- merge(gt1_dt, gt2_dt, by=c("exprds", "entrezID"), suffixes=c("_hicds1", "_hicds2"), all=FALSE)
      stopifnot(!is.na(merge_dt))
      stopifnot(length(unique(merge_dt$exprds)) == 1)
      
      
      tad_rank_corr <- cor(merge_dt$tad_rank_hicds2,merge_dt$tad_rank_hicds1, method=corMeth )
      tad_adjPval_corr <- cor(merge_dt$tad_adjCombPval_hicds2,merge_dt$tad_adjCombPval_hicds1, method=corMeth )
      
      if(exprds == "TCGAluad_norm_luad" & 
         "LG1_40kb" %in% c(hicds1, hicds2) & "ENCSR489OCU_NCI-H460_40kb"  %in% c(hicds1, hicds2)) {
        
        plot_var="tad_rank"
        for(plot_var in c("tad_rank", "tad_adjCombPval")) {
          
          var1 <- merge_dt[,paste0(plot_var, "_hicds1")]
          var2 <- merge_dt[,paste0(plot_var, "_hicds2")]
          
          outFile <- file.path(outFolder, paste0(exprds, "_", hicds1, "_vs_", hicds2, "_", plot_var, ".", plotType))
          do.call(plotType, list(outFile, height=myHeight, width=myWidth))
          densplot(
            x = var1,
            y = var2,
            cex.axis = axisCex,
            cex.lab = labCex,
            # xlab = paste0(hicds1, "\n", plot_var), 
            # ylab = paste0(hicds2, "\n", plot_var), 
            xlab = paste0(gsub("_40kb", "", hicds1)),
            ylab = paste0(gsub("_40kb", "", hicds2)),
            # sub = paste0(hicds2, " vs. ", hicds1),
            main = paste0(exprds)
          )
          mtext(side=3, text = paste0(plot_var), font=3)
          curve(1*x, add=TRUE, lty=2, col="black", lwd=1.5)
          addCorr(
            x = var1, y = var2,
            legPos = "topleft", bty="n"
          )
          foo <- dev.off()
          cat(paste0("... written: ", outFile, "\n"))
          
        }
      }
      data.frame(
        exprds = exprds,
        hicds1 = hicds1,
        hicds2 = hicds2,
        nTADs1 = nrow(gt1_dt),
        nTADs2 = nrow(gt2_dt),
        nIntersectTADs = nrow(merge_dt),
        tad_rank_corr=tad_rank_corr,
        tad_adjPval_corr=tad_adjPval_corr,
        stringsAsFactors = FALSE
      )
    } # end iterating over hicds
    exprds_dt
  } # end iterating over exprdse
  
  
  outFile <- file.path(outFolder, "all_conds_dt.Rdata")
  save(all_conds_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  } else{
    
    outFile <- file.path(outFolder, "all_conds_dt.Rdata")
        
    all_conds_dt <- get(load(outFile))
  }

# all_conds_dt <- get(load("CMP_SAME_CONDITIONS_V2/all_conds_dt.Rdata"))

plot_var="tad_rank"

for(plot_var in c("tad_rank", "tad_adjPval")) {
  
  
  plotTit <- paste0("Same exprds - diff. exprds (n=",nrow(all_conds_dt) , ")")
  mySub <- paste0("corr. of ", plot_var, " (merged on gene ID)")
  
  p2 <- ggdensity(all_conds_dt,
                  fill = "lightgray",
                  x = paste0(plot_var, "_corr"),
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0(corMeth, " corr. for ", plot_var),
                  # add = "median",                  # Add median line. 
                  rug = TRUE,                      # Add marginal rug
                  # color = "signif",
                  # fill = "signif",
                  palette = "jco"
  ) +
    # scale_color_manual(values=my_cols)+
    # scale_fill_manual(values=my_cols)  +
    ggtitle(plotTit, subtitle = mySub)+
    # labs(color=paste0(legTitle),fill=paste0(legTitle)
    labs( y="Density") + 
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
    ) 
  
  outFile <- file.path(outFolder, paste0("sameExprds_diffHicds_", plot_var, "_corr_density.", plotType))
  ggsave(p2, file=outFile, height=myHeight, width=myWidth*1.2)
  cat(paste0("... written: ", outFile, "\n"))
  
}
  
  
  
  
  
  
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




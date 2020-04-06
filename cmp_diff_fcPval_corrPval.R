
# Rscript cmp_diff_fcPval_corrPval.R

require(foreach)
require(doMC)
registerDoMC(40)

# COMP logFC pvals vs corr pval for all CMPs

outFolder <- "CMP_DIFF_FCPVAL_CORRPVAL"
dir.create(outFolder, recursive = TRUE)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myHeight <- myWidth <- 400
plotCex <- 1.4

# v0 = logFC empPval from g2t 100'000 permut; corr empPval from meanCorr across BD permut
# sameNbr = logFC empPval AND corr empPval from meanCorr across BD permut
# random = logFC empPval AND corr empPval from pipeline run for RANDOMMIDPOSSTRICT


all_patterns <- c("", "sameNbr", "random", "randomRescaled", "sameNbrDouble")
all_patterns <- c("sameNbrDouble")
curr_patt = all_patterns[1]

buildTable <- TRUE

pipFolder <- file.path( "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(file.path(pipFolder))
all_hicds1 <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_hicds2 <- all_hicds[grepl("RANDOMMIDPOSSTRICT", all_hicds) | grepl("RANDOMMIDPOSDISC", all_hicds)]
all_hicds <- c(all_hicds1, all_hicds2)
all_hicds <- c(all_hicds1)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

rd_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMSHIFT", "PERMUTG2T", "RANDOMMIDPOSDISC" , "RANDOMMIDPOSSTRICT")

hicds="LG1_40kb"
exprds="TCGAluad_norm_luad"

if(buildTable){
  
  all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      patt_dt <- foreach(curr_patt=all_patterns, .combine='rbind') %do% {
        
        # 9random_runEmpPvalMeanTADLogFC/emp_pval_meanLogFC.Rdata, 9sameNbr_runEmpPvalMeanTADLogFC/emp_pval_meanLogFC.Rdata
        fc_file <- file.path(pipFolder, hicds, exprds,paste0("9", curr_patt, "_runEmpPvalMeanTADLogFC"), "emp_pval_meanLogFC.Rdata")
        if(!file.exists(fc_file)) cat(paste0(fc_file, "\n"))
        fc_pval <- get(load(fc_file))
        # 10random_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata  10random_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata
        if(curr_patt == "" | curr_patt == "sameNbrDouble") {
          corr_pval <- get(load(file.path(pipFolder, hicds, exprds, paste0("10sameNbr", "_runEmpPvalMeanTADCorr"), "emp_pval_meanCorr.Rdata")))
        } else if(curr_patt == "randomRescaled") {
          corr_pval <- get(load(file.path(pipFolder, hicds, exprds, paste0("10random", "_runEmpPvalMeanTADCorr"), "emp_pval_meanCorr.Rdata")))
        } else {
          corr_pval <- get(load(file.path(pipFolder, hicds, exprds, paste0("10", curr_patt, "_runEmpPvalMeanTADCorr"), "emp_pval_meanCorr.Rdata")))
        }# 11random_runEmpPvalCombined/emp_pval_combined.Rdata 11random_runEmpPvalCombined/emp_pval_combined.Rdata
        if(curr_patt == "sameNbr") {
          comb_pval <- get(load(file.path(pipFolder, hicds, exprds, "11sameNbrSameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")))
		}else if(curr_patt == "sameNbrDouble") {
          comb_pval <- get(load(file.path(pipFolder, hicds, exprds, "11sameNbrSameNbrDouble_runEmpPvalCombined", "emp_pval_combined.Rdata")))
		}else if(curr_patt == "") {
          comb_pval <- get(load(file.path(pipFolder, hicds, exprds, "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")))
        }else {
          comb_pval <- get(load(file.path(pipFolder, hicds, exprds,paste0("11", curr_patt, "_runEmpPvalCombined"), "emp_pval_combined.Rdata")))
        }
        
        
        all_regs <- names(fc_pval)
        stopifnot(setequal(names(corr_pval), all_regs))
        stopifnot(setequal(names(comb_pval), all_regs))
        stopifnot(length(fc_pval) == length(corr_pval))
        stopifnot(length(fc_pval) == length(comb_pval))
        
        hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                          gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                               gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                    gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC",
                                         gsub(".+RANDOMMIDPOSSTRICT_40kb", "RANDOMMIDPOSSTRICT",
                                              gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", hicds))))))
        
        if(! hicds_lab %in% rd_patterns) hicds_lab <- "OBSERVED"
        
        
        out_dt <- data.frame(
          hicds=hicds,
          hicds_lab = hicds_lab,
          exprds=exprds,
          region = all_regs,
          pvalType = curr_patt,
          fc_empPval = as.numeric(fc_pval[all_regs]),
          corr_empPval = as.numeric(corr_pval[all_regs]),
          combPval = as.numeric(comb_pval[all_regs]),
          stringsAsFactors = FALSE
        )

        out_dt
        
      }
    patt_dt
      
    }
    exprds_dt
  }
  
  outFile <- file.path(outFolder, paste0("all_result_dt.Rdata"))
  save(all_result_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  
  outFile <- file.path(outFolder, paste0("all_result_dt.Rdata"))
  all_result_dt <- get(load(outFile))
  
}




ds_types <- unique(all_result_dt$hicds_lab)

pval_types <- unique(all_result_dt$pvalType)

for(ds_type in ds_types) {
  
  for(pval_type in pval_types) {
    
    
    plot_dt <- all_result_dt[all_result_dt$hicds_lab == ds_type & all_result_dt$pvalType == pval_type,]
    
    nDS <- length(unique(file.path(plot_dt$hicds, plot_dt$exprds)))
    
    
    
    
    for(i_col in 6:(ncol(plot_dt)-1)) {
      
      
      for(j_col in (i_col+1):ncol(plot_dt)){
        
        my_xlab <- paste0(colnames(plot_dt)[i_col])
        cmp_x <- gsub("adjPvalComb_", "", my_xlab)
        my_xlab <- paste0(my_xlab, " [-log10]")
        
        my_ylab <- paste0(colnames(plot_dt)[j_col])
        cmp_y <- gsub("adjPvalComb_", "", my_ylab)
        my_ylab <- paste0(my_ylab, " [-log10]")
        
        
        cat("... plotting: ", pval_type, " - ", ds_type, " - ", cmp_x, " - ", cmp_y, "\n")
        
        
        my_x <- -log10(plot_dt[,i_col])
        
        
        my_y <- -log10(plot_dt[,j_col])
        
        
        outFile <- file.path(outFolder, paste0(cmp_y, "_vs_", cmp_x, "_adjPvalComb_log10_", pval_type, "_", ds_type, ".", plotType))
        do.call(plotType, list(outFile, height=myHeight, width=myWidth))
        densplot(
          x=my_x,
          y=my_y,
          xlab=my_xlab,
          ylab=my_ylab,
          main=paste0(cmp_y, " vs. ", cmp_x, " - ", pval_type, " - ", ds_type),
          cex.axis=plotCex,
          cex.lab=plotCex,
          cex.main=1
        )
        mtext(side=3, text=paste0("all DS - n=", nDS, "; # TADs=", nrow(plot_dt)))
        addCorr(x=my_x, y  = my_y, bty="n", legPos = "topleft")
        foo <- dev.off()
        cat(paste0("... written: ", outFile, "\n"))
        
        
        
        
      }
      
      
    }
    
    
    
    
  }
  
  
}




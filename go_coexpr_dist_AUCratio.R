startTime <- Sys.time()
cat(paste0("> Rscript go_coexpr_dist_AUCratio.R\n"))

buildTable <- TRUE

# Rscript go_coexpr_dist_AUCratio.R

options(scipen=100)
require(ggpubr)
require(reshape2)
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4


myHeightGG <- 7
myWidthGG <- myHeightGG*1.2


col1 <- get_palette("Dark2", 4)[1]
col2 <- get_palette("Dark2", 4)[2]
col3 <- get_palette("Dark2", 4)[3]
col4 <- get_palette("Dark2", 4)[4]





hicds="Panc1_rep12_40kb"
exprds="TCGApaad_wt_mutKRAS"

outFolder <- file.path("GO_COEXPR_DIST_AUCRATIO")
dir.create(outFolder, recursive = TRUE)

setDir="/media/electron"
setDir=""

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))

pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))

aucFolder <- file.path(mainFolder, "COEXPR_DIST_SAME_GO_SORTNODUP")
stopifnot(dir.exists(aucFolder))

all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds




  # all_hicds=all_hicds[1]
  hicds = all_hicds[1]
  all_auc_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    
    
    # all_exprds = all_exprds[[paste0(hicds)]][1]
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat(paste0("... start ", hicds, " - ", exprds, "\n"))
      
      
      aucFile  <- file.path(aucFolder, hicds, exprds, paste0(hicds,"_", exprds, "_all_auc.Rdata") )
      stopifnot(file.exists(aucFile))
      all_auc <- get(load(aucFile))
      out_dt <- data.frame(all_auc)
      all_cols <- colnames(out_dt)
      out_dt$hicds <- hicds
      out_dt$exprds <- exprds
      out_dt[,c("hicds", "exprds", all_cols)]
      
      
    }# end-foreach iterating over exprds
    
    exprds_dt
    
  }# end-foreach iterating over hicds
  

all_auc_dt$sameTAD_diffTAD_aucRatio <- all_auc_dt$auc_sameTAD_distVect/  all_auc_dt$auc_diffTAD_distVect
all_auc_dt$sameGO_sameTAD_diffTAD_aucRatio <- all_auc_dt$auc_sameGO_sameTAD_distVect/  all_auc_dt$auc_sameGO_diffTAD_distVect



plot_dt <- all_auc_dt[,c("hicds", "exprds", "sameTAD_diffTAD_aucRatio", "sameGO_sameTAD_diffTAD_aucRatio")]
plot_dt_m <- melt(plot_dt, id=c("hicds", "exprds"))

p <- ggdensity(plot_dt_m, 
               title = paste0("coexpr ~ dist AUC ratio (sameTAD/diffTAD)"),
               subtitle=paste0(""),
               x = "value", 
               color = "variable", fill = "variable",
               # add = "mean", rug = TRUE,
               xlab = "AUC ratio",
               palette = c(col1, col2))
p <- p + geom_vline(xintercept = 1, linetype=2)

outFile <- file.path(outFolder, paste0("all_ds_sameGO_AUC_ratio", "_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))





plot_dt <- all_auc_dt[,c("hicds", "exprds", "auc_diffTAD_distVect", "auc_sameTAD_distVect","auc_sameGO_diffTAD_distVect", "auc_sameGO_sameTAD_distVect" )]
plot_dt_m <- melt(plot_dt, id=c("hicds", "exprds"))

p <- ggdensity(plot_dt_m, 
               title = paste0("coexpr ~ dist AUC"),
               subtitle=paste0(""),
               x = "value", 
               color = "variable", fill = "variable",
               # add = "mean", rug = TRUE,
               xlab = "AUC",
               palette = c(col1, col2, col3, col4))

outFile <- file.path(outFolder, paste0("all_ds_sameGO_AUC", "_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))



######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


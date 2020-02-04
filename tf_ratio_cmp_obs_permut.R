startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)
require(dplyr)

# Rscript tf_ratio_cmp_obs_permut.R crisp
# Rscript tf_ratio_cmp_obs_permut.R chea3_all
# Rscript tf_ratio_cmp_obs_permut.R chea3_lung 
# Rscript tf_ratio_cmp_obs_permut.R trrust
# Rscript tf_ratio_cmp_obs_permut.R motifmap
# 
#

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- ifelse(plotType=="png", 600, 9)
plotCex <- 1.4


nPermut <- 1000

tadSignifThresh <- 0.01

dsIn <- "chea3_lung"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1 | length(args) == 3)
dsIn <- args[1]
if(length(args) == 3) {
  all_hicds <- args[2]
  all_exprds <- args[3]
  
} else {
  all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))
}


stopifnot(dsIn %in% c("crisp",  "trrust", "motifmap", "chea3_all", "chea3_lung"))

outFolder <- file.path(paste0("TF_RATIO_CMP_OBS_PERMUT_", toupper(dsIn)))
dir.create(outFolder, recursive = TRUE)


final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
final_dt$id <- paste0(final_dt$hicds, ".", final_dt$exprds, ".", final_dt$region)
signif_id <- final_dt$id[final_dt$adjPvalComb <= tadSignifThresh]

#######################################
############# OBSERVED DATA
#######################################

all_obs_data <- get(load(file.path(paste0("TF_RATIO_OBS_", toupper(dsIn)), "all_ratio_data.Rdata")))

all_obs_nGenes <- unlist(lapply(all_obs_data, function(sl) lapply(sl[[2]], function(x) x[[paste0("exprds_nGenes")]])))
all_obs_nTargetByGenes <- unlist(lapply(all_obs_data, function(sl) lapply(sl[[2]], function(x) x[[paste0("exprds_nTarget")]]/ x[["exprds_nGenes"]] )))
all_obs_nTargetExprByGenes <- unlist(lapply(all_obs_data, function(sl) lapply(sl[[2]], function(x) x[[paste0("exprds_nTargetExpr")]]/ x[["exprds_nGenes"]] )))
all_obs_nTargetWithBDByGenes <- unlist(lapply(all_obs_data, function(sl) lapply(sl[[2]], function(x) x[[paste0("exprds_nTargetWithBD")]]/ x[["exprds_nGenes"]] )))
all_obs_nTargetWithBDExprByGenes <- unlist(lapply(all_obs_data, function(sl) lapply(sl[[2]], function(x) x[[paste0("exprds_nTargetWithBDExpr")]]/ x[["exprds_nGenes"]] )))
all_obs_nRegByGenes <-unlist(lapply(all_obs_data, function(sl) lapply(sl[[2]], function(x) x[[paste0("exprds_nReg")]]/ x[["exprds_nGenes"]] )))

stopifnot(length(all_obs_nGenes) == length(all_obs_nTargetByGenes))
stopifnot(length(all_obs_nGenes) == length(all_obs_nTargetExprByGenes))
stopifnot(length(all_obs_nGenes) == length(all_obs_nTargetWithBDByGenes))
stopifnot(length(all_obs_nGenes) == length(all_obs_nTargetWithBDExprByGenes))
stopifnot(length(all_obs_nGenes) == length(all_obs_nRegByGenes))

stopifnot(names(all_obs_nGenes) %in% final_dt$id)

all_obsSignif_nGenes <- all_obs_nGenes[names(all_obs_nGenes) %in% signif_id]
all_obsSignif_nTargetByGenes <- all_obs_nTargetByGenes[names(all_obs_nTargetByGenes) %in% signif_id] 
all_obsSignif_nTargetExprByGenes <- all_obs_nTargetExprByGenes[names(all_obs_nTargetExprByGenes) %in% signif_id]
all_obsSignif_nTargetWithBDByGenes <- all_obs_nTargetWithBDByGenes[names(all_obs_nTargetWithBDByGenes) %in% signif_id]
all_obsSignif_nTargetWithBDExprByGenes <- all_obs_nTargetWithBDExprByGenes[names(all_obs_nTargetWithBDExprByGenes) %in% signif_id]
all_obsSignif_nRegByGenes <- all_obs_nRegByGenes[names(all_obs_nRegByGenes) %in% signif_id]

 
all_obsNotSignif_nGenes <- all_obs_nGenes[!names(all_obs_nGenes) %in% signif_id]
all_obsNotSignif_nTargetByGenes <- all_obs_nTargetByGenes[!names(all_obs_nTargetByGenes) %in% signif_id] 
all_obsNotSignif_nTargetExprByGenes <- all_obs_nTargetExprByGenes[!names(all_obs_nTargetExprByGenes) %in% signif_id]
all_obsNotSignif_nTargetWithBDByGenes <- all_obs_nTargetWithBDByGenes[!names(all_obs_nTargetWithBDByGenes) %in% signif_id]
all_obsNotSignif_nTargetWithBDExprByGenes <- all_obs_nTargetWithBDExprByGenes[!names(all_obs_nTargetWithBDExprByGenes) %in% signif_id]
all_obsNotSignif_nRegByGenes <- all_obs_nRegByGenes[!names(all_obs_nRegByGenes) %in% signif_id]


#######################################
############# PREP DATA PERMUT-CORR
#######################################

all_permutCorr_data <- get(load(file.path(paste0("TF_RATIO_CORRPERMUT_", toupper(dsIn)), "all_sample_ratio_data.Rdata")))

all_permutCorr_nGenes <- unlist(lapply(all_permutCorr_data, function(sl) lapply(sl,  function(x) x[[paste0("exprds_nGenes")]])))
all_permutCorr_nTargetByGenes <- unlist(lapply(all_permutCorr_data, function(sl) lapply(sl, function(x) x[[paste0("exprds_nTarget")]]/ x[["exprds_nGenes"]] )))
all_permutCorr_nTargetExprByGenes <- unlist(lapply(all_permutCorr_data, function(sl) lapply(sl, function(x) x[[paste0("exprds_nTargetExpr")]]/ x[["exprds_nGenes"]] )))
all_permutCorr_nTargetWithBDByGenes <- unlist(lapply(all_permutCorr_data, function(sl) lapply(sl, function(x) x[[paste0("exprds_nTargetWithBD")]]/ x[["exprds_nGenes"]] )))
all_permutCorr_nTargetWithBDExprByGenes <- unlist(lapply(all_permutCorr_data, function(sl) lapply(sl, function(x) x[[paste0("exprds_nTargetWithBDExpr")]]/ x[["exprds_nGenes"]] )))
all_permutCorr_nRegByGenes <-unlist(lapply(all_permutCorr_data, function(sl) lapply(sl, function(x) x[[paste0("exprds_nReg")]]/ x[["exprds_nGenes"]] )))

stopifnot(length(all_permutCorr_nGenes) == length(all_permutCorr_nTargetByGenes))
stopifnot(length(all_permutCorr_nGenes) == length(all_permutCorr_nTargetExprByGenes))
stopifnot(length(all_permutCorr_nGenes) == length(all_permutCorr_nTargetWithBDByGenes))
stopifnot(length(all_permutCorr_nGenes) == length(all_permutCorr_nTargetWithBDExprByGenes))
stopifnot(length(all_permutCorr_nGenes) == length(all_permutCorr_nRegByGenes))


#######################################
############# COLLECT ALL PERMUT-g2t
#######################################

# all_permutG2t_files=all_permutG2t_files[1]
all_permutG2t_files <- list.files(file.path(paste0("TF_RATIO_G2TPERMUT_", toupper(dsIn))), recursive = TRUE, pattern=paste0("_nPermut", nPermut, "_tf_ratios.Rdata"), full.names = TRUE)


n_hicds <-  gsub("(.+)_TCGA.+", "\\1", basename(all_permutG2t_files)) 
n_exprds <-  gsub(".+_(TCGA.+)_nPermut.+", "\\1", basename(all_permutG2t_files)) 

all_permutG2t_data <- foreach(g2tpermut_file = all_permutG2t_files) %dopar% {
 get(load(g2tpermut_file))
}
names(all_permutG2t_data) <- file.path(n_hicds, n_exprds)



all_permutG2t_nGenes <- unlist(lapply(all_permutG2t_data,  function(permut) lapply(permut, function(x) x[[paste0("exprds_nGenes")]])))
all_permutG2t_nTargetByGenes <- unlist(lapply(all_permutG2t_data, function(permut) lapply(permut, function(x) x[[paste0("exprds_nTarget")]]/ x[["exprds_nGenes"]] )))
all_permutG2t_nTargetExprByGenes <- unlist(lapply(all_permutG2t_data, function(permut) lapply(permut, function(x) x[[paste0("exprds_nTargetExpr")]]/ x[["exprds_nGenes"]] )))
all_permutG2t_nTargetWithBDByGenes <- unlist(lapply(all_permutG2t_data,  function(permut) lapply(permut, function(x) x[[paste0("exprds_nTargetWithBD")]]/ x[["exprds_nGenes"]] )))
all_permutG2t_nTargetWithBDExprByGenes <- unlist(lapply(all_permutG2t_data, function(permut) lapply(permut,  function(x) x[[paste0("exprds_nTargetWithBDExpr")]]/ x[["exprds_nGenes"]] )))
all_permutG2t_nRegByGenes <-unlist(lapply(all_permutG2t_data,  function(permut) lapply(permut, function(x) x[[paste0("exprds_nReg")]]/ x[["exprds_nGenes"]] )))

stopifnot(length(all_permutG2t_nGenes) == length(all_permutG2t_nTargetByGenes))
stopifnot(length(all_permutG2t_nGenes) == length(all_permutG2t_nTargetExprByGenes))
stopifnot(length(all_permutG2t_nGenes) == length(all_permutG2t_nTargetWithBDByGenes))
stopifnot(length(all_permutG2t_nGenes) == length(all_permutG2t_nTargetWithBDExprByGenes))
stopifnot(length(all_permutG2t_nGenes) == length(all_permutG2t_nRegByGenes))




##################################################################

all_vars <- c("nGenes", "nTargetByGenes", "nTargetExprByGenes", "nTargetWithBDByGenes", "nTargetWithBDExprByGenes", "nRegByGenes")
curr_var = all_vars[6]
curr_var = all_vars[1]
for(curr_var in all_vars) {

  cat(paste0("... all ds -  ", curr_var, "\n"))

  all_plot_list <- list(
    get(paste0("all_obs_", curr_var)),
    get(paste0("all_obsSignif_", curr_var)),
    get(paste0("all_obsNotSignif_", curr_var)),
    get(paste0("all_permutCorr_", curr_var)),
    get(paste0("all_permutG2t_", curr_var))
  )

  names(all_plot_list) <- c(paste0("all_obs_", curr_var),
                        paste0("all_obsSignif_", curr_var), 
                        paste0("all_obsNotSignif_", curr_var), 
                        paste0("all_permutCorr_", curr_var), paste0("all_permutG2t_", curr_var))

  save(all_plot_list, file="all_plot_list.Rdata", version=2)
  
  plot_list <- Filter(all_plot_list, f = function(x) length(x)>0 )

  outFile <- file.path(outFolder, paste0("allDS_", curr_var, "_multiDens.", plotType ))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(
    plot_list,
    plotTit = paste0(curr_var)
  )
  mtext(side=3, text = paste0("all data"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))


}


#***************************************************************************************************************************
### DATASET-SPECIFIC
#***************************************************************************************************************************


for(i_ds in 1:length(n_hicds)) {
  
  hicds <- n_hicds[i_ds]
  exprds <- n_exprds[i_ds]
  
  cat(paste0(hicds, " - ", exprds, "\n"))
  
  # ds_obs_nGenes <- unlist(lapply(all_obs_data[[paste0(hicds)]][[2]], function(x) x[["exprds_nGenes"]] ))
  # ds_obs_nTargetByGenes <- unlist(lapply(all_obs_data[[paste0(hicds)]][[2]], function(x) x[[paste0("exprds_nTarget")]]/ x[["exprds_nGenes"]] ))
  # ds_obs_nTargetExprByGenes <- unlist(lapply(all_obs_data[[paste0(hicds)]][[2]], function(x) x[[paste0("exprds_nTargetExpr")]]/ x[["exprds_nGenes"]] ))
  # ds_obs_nTargetWithBDByGenes <- unlist(lapply(all_obs_data[[paste0(hicds)]][[2]],  function(x) x[[paste0("exprds_nTargetWithBD")]]/ x[["exprds_nGenes"]] ))
  # ds_obs_nTargetWithBDExprByGenes <- unlist(lapply(all_obs_data[[paste0(hicds)]][[2]], function(x) x[[paste0("exprds_nTargetWithBDExpr")]]/ x[["exprds_nGenes"]] ))
  # ds_obs_nRegByGenes <- unlist(lapply(all_obs_data[[paste0(hicds)]][[2]],  function(x) x[[paste0("exprds_nReg")]]/ x[["exprds_nGenes"]] ))
  
  ds_obs_nGenes <- unlist(all_obs_data[[paste0(hicds)]][[2]][[paste0(exprds)]][["exprds_nGenes"]] )
  ds_obs_nTargetByGenes <- unlist(all_obs_data[[paste0(hicds)]][[2]][[paste0(exprds)]][[paste0("exprds_nTarget")]]/ ds_obs_nGenes )
  ds_obs_nTargetExprByGenes <- unlist(all_obs_data[[paste0(hicds)]][[2]][[paste0(exprds)]][[paste0("exprds_nTargetExpr")]]/ ds_obs_nGenes )
  ds_obs_nTargetWithBDByGenes <- unlist(all_obs_data[[paste0(hicds)]][[2]][[paste0(exprds)]][[paste0("exprds_nTargetWithBD")]]/ ds_obs_nGenes)
  ds_obs_nTargetWithBDExprByGenes <- unlist(all_obs_data[[paste0(hicds)]][[2]][[paste0(exprds)]][[paste0("exprds_nTargetWithBDExpr")]]/ds_obs_nGenes )
  ds_obs_nRegByGenes <- unlist(all_obs_data[[paste0(hicds)]][[2]][[paste0(exprds)]][[paste0("exprds_nReg")]]/ ds_obs_nGenes )
  
  stopifnot(length(ds_obs_nGenes) == length(ds_obs_nTargetByGenes))
  stopifnot(length(ds_obs_nGenes) == length(ds_obs_nTargetExprByGenes))
  stopifnot(length(ds_obs_nGenes) == length(ds_obs_nTargetWithBDByGenes))
  stopifnot(length(ds_obs_nGenes) == length(ds_obs_nTargetWithBDExprByGenes))
  stopifnot(length(ds_obs_nGenes) == length(ds_obs_nRegByGenes))
  
  
  ds_obsSignif_nGenes <- ds_obs_nGenes[names(ds_obs_nGenes) %in% signif_id]
  ds_obsSignif_nTargetByGenes <- ds_obs_nTargetByGenes[names(ds_obs_nTargetByGenes) %in% signif_id] 
  ds_obsSignif_nTargetExprByGenes <- ds_obs_nTargetExprByGenes[names(ds_obs_nTargetExprByGenes) %in% signif_id]
  ds_obsSignif_nTargetWithBDByGenes <- ds_obs_nTargetWithBDByGenes[names(ds_obs_nTargetWithBDByGenes) %in% signif_id]
  ds_obsSignif_nTargetWithBDExprByGenes <- ds_obs_nTargetWithBDExprByGenes[names(ds_obs_nTargetWithBDExprByGenes) %in% signif_id]
  ds_obsSignif_nRegByGenes <- ds_obs_nRegByGenes[names(ds_obs_nRegByGenes) %in% signif_id]
  
  
  ds_obsNotSignif_nGenes <- ds_obs_nGenes[!names(ds_obs_nGenes) %in% signif_id]
  ds_obsNotSignif_nTargetByGenes <- ds_obs_nTargetByGenes[!names(ds_obs_nTargetByGenes) %in% signif_id] 
  ds_obsNotSignif_nTargetExprByGenes <- ds_obs_nTargetExprByGenes[!names(ds_obs_nTargetExprByGenes) %in% signif_id]
  ds_obsNotSignif_nTargetWithBDByGenes <- ds_obs_nTargetWithBDByGenes[!names(ds_obs_nTargetWithBDByGenes) %in% signif_id]
  ds_obsNotSignif_nTargetWithBDExprByGenes <- ds_obs_nTargetWithBDExprByGenes[!names(ds_obs_nTargetWithBDExprByGenes) %in% signif_id]
  ds_obsNotSignif_nRegByGenes <- ds_obs_nRegByGenes[!names(ds_obs_nRegByGenes) %in% signif_id]
  
  # ds_permutCorr_nGenes <- unlist(lapply(all_permutCorr_data[[paste0(hicds)]], function(x) x[[paste0("exprds_nGenes")]]))
  # ds_permutCorr_nTargetByGenes <- unlist(lapply(all_permutCorr_data[[paste0(hicds)]], function(x) x[[paste0("exprds_nTarget")]]/ x[["exprds_nGenes"]] ))
  # ds_permutCorr_nTargetExprByGenes <- unlist(lapply(all_permutCorr_data[[paste0(hicds)]], function(x) x[[paste0("exprds_nTargetExpr")]]/ x[["exprds_nGenes"]] ))
  # ds_permutCorr_nTargetWithBDByGenes <- unlist(lapply(all_permutCorr_data[[paste0(hicds)]],  function(x) x[[paste0("exprds_nTargetWithBD")]]/ x[["exprds_nGenes"]] ))
  # ds_permutCorr_nTargetWithBDExprByGenes <- unlist(lapply(all_permutCorr_data[[paste0(hicds)]], function(x) x[[paste0("exprds_nTargetWithBDExpr")]]/ x[["exprds_nGenes"]] ))
  # ds_permutCorr_nRegByGenes <- unlist(lapply(all_permutCorr_data[[paste0(hicds)]], function(x) x[[paste0("exprds_nReg")]]/ x[["exprds_nGenes"]] ))
  
  
  ds_permutCorr_nGenes <- unlist(all_permutCorr_data[[paste0(hicds)]][[paste0(exprds)]][[paste0("exprds_nGenes")]])
  ds_permutCorr_nTargetByGenes <- unlist(all_permutCorr_data[[paste0(hicds)]][[paste0(exprds)]][[paste0("exprds_nTarget")]]/ ds_permutCorr_nGenes )
  ds_permutCorr_nTargetExprByGenes <- unlist(all_permutCorr_data[[paste0(hicds)]][[paste0(exprds)]][[paste0("exprds_nTargetExpr")]]/ ds_permutCorr_nGenes )
  ds_permutCorr_nTargetWithBDByGenes <- unlist(all_permutCorr_data[[paste0(hicds)]][[paste0(exprds)]][[paste0("exprds_nTargetWithBD")]]/ ds_permutCorr_nGenes )
  ds_permutCorr_nTargetWithBDExprByGenes <- unlist(all_permutCorr_data[[paste0(hicds)]][[paste0(exprds)]][[paste0("exprds_nTargetWithBDExpr")]]/ ds_permutCorr_nGenes)
  ds_permutCorr_nRegByGenes <- unlist(all_permutCorr_data[[paste0(hicds)]][[paste0(exprds)]][[paste0("exprds_nReg")]]/ ds_permutCorr_nGenes )
  
  stopifnot(length(ds_permutCorr_nGenes) == length(ds_permutCorr_nTargetByGenes))
  stopifnot(length(ds_permutCorr_nGenes) == length(ds_permutCorr_nTargetExprByGenes))
  stopifnot(length(ds_permutCorr_nGenes) == length(ds_permutCorr_nTargetWithBDByGenes))
  stopifnot(length(ds_permutCorr_nGenes) == length(ds_permutCorr_nTargetWithBDExprByGenes))
  stopifnot(length(ds_permutCorr_nGenes) == length(ds_permutCorr_nRegByGenes))
  
  
  ds_g2tPermut <- all_permutG2t_data[[file.path(hicds, exprds)]]
  ds_permutG2t_nGenes <- unlist(lapply(ds_g2tPermut,  function(x) x[[paste0("exprds_nGenes")]]))
  ds_permutG2t_nTargetByGenes <- unlist(lapply(ds_g2tPermut,  function(x) x[[paste0("exprds_nTarget")]]/ x[["exprds_nGenes"]] ))
  ds_permutG2t_nTargetExprByGenes <- unlist(lapply(ds_g2tPermut,  function(x) x[[paste0("exprds_nTargetExpr")]]/ x[["exprds_nGenes"]] ))
  ds_permutG2t_nTargetWithBDByGenes <- unlist(lapply(ds_g2tPermut,   function(x) x[[paste0("exprds_nTargetWithBD")]]/ x[["exprds_nGenes"]] ))
  ds_permutG2t_nTargetWithBDExprByGenes <- unlist(lapply(ds_g2tPermut,  function(x) x[[paste0("exprds_nTargetWithBDExpr")]]/ x[["exprds_nGenes"]] ))
  ds_permutG2t_nRegByGenes <-unlist(lapply(ds_g2tPermut,function(x) x[[paste0("exprds_nReg")]]/ x[["exprds_nGenes"]] ))
  
  stopifnot(length(ds_permutG2t_nGenes) == length(ds_permutG2t_nTargetByGenes))
  stopifnot(length(ds_permutG2t_nGenes) == length(ds_permutG2t_nTargetExprByGenes))
  stopifnot(length(ds_permutG2t_nGenes) == length(ds_permutG2t_nTargetWithBDByGenes))
  stopifnot(length(ds_permutG2t_nGenes) == length(ds_permutG2t_nTargetWithBDExprByGenes))
  stopifnot(length(ds_permutG2t_nGenes) == length(ds_permutG2t_nRegByGenes))
  
  
  all_vars <- c("nGenes", "nTargetByGenes", "nTargetExprByGenes", "nTargetWithBDByGenes", "nTargetWithBDExprByGenes", "nRegByGenes")
  curr_var = all_vars[6]
  curr_var = all_vars[1]
  for(curr_var in all_vars) {
    
    cat(paste0("... ", exprds,  " - ", hicds, " - ", curr_var, "\n"))
    
    plot_list_all <- list(
      get(paste0("ds_obs_", curr_var)),
      get(paste0("ds_obsSignif_", curr_var)),
      get(paste0("ds_obsNotSignif_", curr_var)),
      get(paste0("ds_permutCorr_", curr_var)),
      get(paste0("ds_permutG2t_", curr_var))
    )
    names(plot_list_all) <- c(paste0("ds_obs_", curr_var), 
                          paste0("ds_obsSignif_", curr_var),
                          paste0("ds_obsNotSignif_", curr_var),
                          paste0("ds_permutCorr_", curr_var), paste0("ds_permutG2t_", curr_var))
    
    save(plot_list_all, file="plot_list_all.Rdata", version=2)
    
    plot_list <- Filter(plot_list_all, f = function(x) length(x)>0 )
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", curr_var, "_multiDens.", plotType ))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens(
      plot_list,
      plotTit = paste0(curr_var)
    )
    mtext(side=3, text = paste0(hicds, " - ", exprds))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
  }
  
  
}






#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))

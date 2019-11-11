startTime <- Sys.time()
cat(paste0("> Rscript paralogs_and_tads_obsOnly.R\n"))

script_name <- "paralogs_and_tads_obsOnly.R"


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript paralogs_and_tads_obsOnly.R
buildTable <- FALSE

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

source("paralog_utils_fct.R")
nRandom <- 1000
# nRandom <- 5
stopifnot(nRandom>1)
set.seed(10112019)

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(pipFolder)
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

outFolder <- file.path("PARALOGS_AND_TADS_OBSONLY")
dir.create(outFolder, recursive = TRUE)

################################################################################################################################################
######################################################### PREP PARALOG DATA
################################################################################################################################################

id_dt <- read.delim("ensembl_entrez_biomart.txt", sep=",", stringsAsFactors=FALSE, header=TRUE)

paralog_dt <- read.delim("ensembl_paralogs_biomart.txt", sep=",", stringsAsFactors=FALSE, header=TRUE)
paralog_dt <- paralog_dt[!paralog_dt$Human.paralogue.gene.stable.ID=="",]

par_dt_sub <- paralog_dt[, c("Gene.stable.ID", "Human.paralogue.gene.stable.ID" )]
id_dt_sub <- id_dt[,c("Gene.stable.ID", "EntrezGene.ID")]
id_dt_sub$EntrezGene.ID <- as.character(id_dt_sub$EntrezGene.ID)

out_dt_tmp <- merge(par_dt_sub, id_dt_sub, by="Gene.stable.ID")

colnames(id_dt_sub)[colnames(id_dt_sub) == "Gene.stable.ID"] <- "Human.paralogue.gene.stable.ID"
colnames(id_dt_sub)[colnames(id_dt_sub) == "EntrezGene.ID"] <- "Human.paralogue.entrezGene.ID"

out_dt <- merge(out_dt_tmp, id_dt_sub, by="Human.paralogue.gene.stable.ID")
entrez_paraDT <- out_dt[,c("EntrezGene.ID", "Human.paralogue.entrezGene.ID")]
entrez_paraDT <- na.omit(entrez_paraDT)
nrow(entrez_paraDT)
# 200196
################################################################################################################################################


####################################################################################################################################### >>> prepare the data
if(buildTable) {
  cat("... start buildTable \n")
  
  hicds = all_hicds[1]
  
  # all_hicds = all_hicds[1:2]
  
  all_ds_results <- foreach(hicds = all_hicds) %dopar% {
    cat(paste0("... start ", hicds, "\n"))  
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    tad_file <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt") #  used for permutation
    stopifnot(file.exists(hicds_file))
    stopifnot(file.exists(tad_file))
    obs_g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    obs_g2t_dt$entrezID <- as.character(obs_g2t_dt$entrezID)
    stopifnot(nrow(obs_g2t_dt) > 0 )
    
    av_entrez <- unique(c(entrez_paraDT$EntrezGene.ID[entrez_paraDT$EntrezGene.ID %in% obs_g2t_dt$entrezID & 
                                                        entrez_paraDT$Human.paralogue.entrezGene.ID %in% obs_g2t_dt$entrezID],
                    entrez_paraDT$Human.paralogue.entrezGene.ID[entrez_paraDT$EntrezGene.ID %in% obs_g2t_dt$entrezID &
                                                  entrez_paraDT$Human.paralogue.entrezGene.ID %in% obs_g2t_dt$entrezID]))
    
    
              para_g2t_dt_tmp <- inner_join(entrez_paraDT, obs_g2t_dt[,c("entrezID", "region")], by = c("EntrezGene.ID" = "entrezID"))
              colnames(para_g2t_dt_tmp)[colnames(para_g2t_dt_tmp) == "region"] <- "ref_region"
              para_g2t_dt <- inner_join(para_g2t_dt_tmp, obs_g2t_dt[,c("entrezID", "region")], by = c("Human.paralogue.entrezGene.ID" = "entrezID"))
              colnames(para_g2t_dt)[colnames(para_g2t_dt) == "region"] <- "paralog_region"
              nrow(para_g2t_dt)
              # 158274
              stopifnot(is.character(para_g2t_dt$EntrezGene.ID))
              stopifnot(is.character(para_g2t_dt$Human.paralogue.entrezGene.ID))
              stopifnot(setequal(av_entrez, unique(c(para_g2t_dt$EntrezGene.ID, para_g2t_dt$Human.paralogue.entrezGene.ID))))
              
              obs_sameTAD_ratio <- sum(para_g2t_dt$ref_region == para_g2t_dt$paralog_region)/nrow(para_g2t_dt)
              
              cat(paste0("... create paralog sets:\t", Sys.time() , " - "))
              para_list_2 <- by(para_g2t_dt, para_g2t_dt$EntrezGene.ID, function(sub_dt) sort(unique(c(sub_dt$EntrezGene.ID, sub_dt$Human.paralogue.entrezGene.ID)) ))
              cat(paste0(Sys.time() , " \n "))
              para_list <- para_list_2
              names(para_list) <- NULL
              para_list <- unique(para_list)
              length(para_list) # 3491
              length(para_list_2) # 12993
              
              # > para_list[[35]]
              # [1] "100131539" "100131980" "100132396" "169270"    "728957"
              # > para_list[[39]]
              # [1] "100131539" "100131980" "100132396" "169270"    "440077"    "728957"   
          
              cat(paste0("... remove nested:\t", Sys.time(), " - "))
              
              n_match_paras <- unlist(lapply(para_list, function(mylist) sum(unlist(lapply(para_list, function(x) all(mylist %in% x))))))
              stopifnot(length(n_match_paras) == length(para_list))
              stopifnot(n_match_paras > 0)
              cat(paste0(Sys.time(), "\n"))
              
              para_list_nosub <- para_list[which(n_match_paras == 1)]
              length(para_list_nosub)
          
              cat(paste0("... still duplicated entrez:\t", sum(duplicated(unlist(para_list_nosub))), "\n"))
              
              para_list_tads <- lapply(para_list_nosub, function(sublist) {
                unique(unlist(sapply(sublist, function(sublist_gene) obs_g2t_dt$region[obs_g2t_dt$entrezID == sublist_gene])))
              })
              
              obs_para_list_nGenes <- lengths(para_list_nosub)
              obs_para_list_nTads <- lengths(para_list_tads)
              
              stopifnot(length(obs_para_list_nGenes) == length(obs_para_list_nTads))
              stopifnot(obs_para_list_nGenes >= obs_para_list_nTads)
              
              obs_para_list_ratioTadsGenes <- obs_para_list_nTads/obs_para_list_nGenes
              
              hist(obs_para_list_nTads, breaks=1:max(obs_para_list_nTads))
    
    
    check_fct <- get_paralog_concord_with_tad(dt_para=entrez_paraDT,
                                             dt_g2t=obs_g2t_dt, 
                                             check_entrez=av_entrez)
    
    save(obs_sameTAD_ratio, file=file.path(outFolder, "obs_sameTAD_ratio.Rdata"), version=2)
    save(obs_para_list_nTads, file=file.path(outFolder, "obs_para_list_nTads.Rdata"), version=2)
    save(check_fct, file=file.path(outFolder, "check_fct.Rdata"), version=2)
      
    cat(paste0("check_fct[[sameTAD_ratio]] \t=\t", round(check_fct[["sameTAD_ratio"]], 4), "\n"))
    
    stopifnot(round(obs_sameTAD_ratio,4) == round(check_fct[["sameTAD_ratio"]], 4))
    
    stopifnot(obs_para_list_nTads == check_fct[["paralogs_list_nTads"]])
    stopifnot(obs_para_list_nGenes == check_fct[["paralogs_list_nGenes"]])
    stopifnot(obs_para_list_ratioTadsGenes == check_fct[["paralogs_list_ratioTadsGenes"]])
    
    cat(paste0("stopifnot(obs_para_list_nTads == check_fct[[paralogs_list_nTads]])   \t=\t", 
               stopifnot(obs_para_list_nTads == check_fct[["paralogs_list_nTads"]]), "\n"))  
    
    
    
    list(obs_sameTAD_ratio=obs_sameTAD_ratio,
         obs_para_list_nTads=obs_para_list_nTads,
         obs_para_list_ratioTadsGenes = obs_para_list_ratioTadsGenes
         )
    
    
  } # end-foreach hicds
  names(all_ds_results) <- all_hicds
  outFile <- file.path(outFolder, "all_ds_results.Rdata")
  save(all_ds_results, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
} else {# END-IF BUILDTABLE
  
  outFile <- file.path(outFolder, "all_ds_results.Rdata")
  all_ds_results <- get(load(outFile))
  
}
all_obs_ratio <- unlist(lapply(all_ds_results, function(x) x[["obs_sameTAD_ratio"]]))
# all_rd_ratio <- unlist(lapply(all_ds_results, function(sublist) lapply(sublist[["permut_results"]], function(x) x[["sameTAD_ratio"]])))

outFile <- file.path(outFolder, paste0("ratioParalogPairsInSameTAD.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth*1.2))  
plot_multiDens(
  list(
    obs_ratio = all_obs_ratio
    # permut_ratio = all_rd_ratio
  ),
  plotTit="ratio of paralog pairs in same TAD"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

all_obs_ratioTadsGenes <- unlist(lapply(all_ds_results, function(x)x[["obs_para_list_ratioTadsGenes"]]))
# all_rd_nTads <- unlist(lapply(all_ds_results, function(sublist) lapply(sublist[["permut_results"]], function(x) x[["paralogs_list_nTads"]])))


outFile <- file.path(outFolder, paste0("ratioTadsGenes_paralogGroup.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth*1.2))  
plot_multiDens(
  list(
    obs_nTads = all_obs_ratioTadsGenes
    # permut_nTads = all_rd_nTads
  ),
  plotTit="ratio of #TADs/#genes paralog groups"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


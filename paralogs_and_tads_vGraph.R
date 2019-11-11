startTime <- Sys.time()
cat(paste0("> Rscript paralogs_and_tads_vGraph.R\n"))

script_name <- "paralogs_and_tads_vGraph.R"

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
require(igraph)

# Rscript paralogs_and_tads_vGraph.R
buildTable <- TRUE

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

source("paralog_utils_fct_vGraph.R")
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

outFolder <- file.path("PARALOGS_AND_TADS_vGraph_all")
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
  
  # all_hicds = all_hicds[1]
  
  all_ds_results <- foreach(hicds = all_hicds) %do% {
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
    
    
              # para_g2t_dt_tmp <- inner_join(entrez_paraDT, obs_g2t_dt[,c("entrezID", "region")], by = c("EntrezGene.ID" = "entrezID"))
              # colnames(para_g2t_dt_tmp)[colnames(para_g2t_dt_tmp) == "region"] <- "ref_region"
              # para_g2t_dt <- inner_join(para_g2t_dt_tmp, obs_g2t_dt[,c("entrezID", "region")], by = c("Human.paralogue.entrezGene.ID" = "entrezID"))
              # colnames(para_g2t_dt)[colnames(para_g2t_dt) == "region"] <- "paralog_region"
              # nrow(para_g2t_dt)
              # # 158274
              # stopifnot(is.character(para_g2t_dt$EntrezGene.ID))
              # stopifnot(is.character(para_g2t_dt$Human.paralogue.entrezGene.ID))
              # stopifnot(setequal(av_entrez, unique(c(para_g2t_dt$EntrezGene.ID, para_g2t_dt$Human.paralogue.entrezGene.ID))))
              # 
              # obs_sameTAD_ratio <- sum(para_g2t_dt$ref_region == para_g2t_dt$paralog_region)/nrow(para_g2t_dt)
              # 
              # cat(paste0("... create paralog sets:\t", Sys.time() , " - "))
              # 
              # para_g <- graph_from_data_frame(para_g2t_dt[,c("EntrezGene.ID", "Human.paralogue.entrezGene.ID")], directed = FALSE, vertices = NULL)
              # para_list_nosub <- lapply(decompose.graph(para_g), function(x) names(V(x)))
              # 
              # length(para_list_nosub)
              # 
              # cat(paste0("... still duplicated entrez:\t", sum(duplicated(unlist(para_list_nosub))), "\n"))
              # 
              ## para_list_tads <- lapply(para_list_nosub, function(sublist) {
              ##   unique(unlist(sapply(sublist, function(sublist_gene) obs_g2t_dt$region[obs_g2t_dt$entrezID == sublist_gene])))
              ## })
              
              # para_list_tads <- lapply(para_list_nosub, function(sublist) {
              #   unique(unlist(obs_g2t_dt$region[obs_g2t_dt$entrezID %in% sublist]))
              # })
              # 
              # 
              # obs_para_list_nGenes <- lengths(para_list_nosub)
              # obs_para_list_nTads <- lengths(para_list_tads)
              # 
              # stopifnot(length(obs_para_list_nGenes) == length(obs_para_list_nTads))
              # stopifnot(obs_para_list_nGenes >= obs_para_list_nTads)
              # 
              # obs_para_list_ratioTadsGenes <- obs_para_list_nTads/obs_para_list_nGenes
              # 
              # hist(obs_para_list_nTads, breaks=1:max(obs_para_list_nTads))
    
    
    check_fct <- get_paralog_concord_with_tad(dt_para=entrez_paraDT,
                                             dt_g2t=obs_g2t_dt, 
                                             check_entrez=av_entrez)
    
    obs_sameTAD_ratio <- check_fct[["sameTAD_ratio"]]
    obs_para_list_nTads <- check_fct[["paralogs_list_nTads"]]
    obs_para_list_ratioTadsGenes <-  check_fct[["paralogs_list_ratioTadsGenes"]]
    # save(obs_sameTAD_ratio, file=file.path(outFolder, "obs_sameTAD_ratio.Rdata"), version=2)
    # save(obs_para_list_nTads, file=file.path(outFolder, "obs_para_list_nTads.Rdata"), version=2)
    # save(check_fct, file=file.path(outFolder, "check_fct.Rdata"), version=2)
      
    cat(paste0("check_fct[[sameTAD_ratio]] \t=\t", round(check_fct[["sameTAD_ratio"]], 4), "\n"))
    
    # stopifnot(round(obs_sameTAD_ratio,4) == round(check_fct[["sameTAD_ratio"]], 4))
    
    # stopifnot(obs_para_list_nTads == check_fct[["paralogs_list_nTads"]])
    # stopifnot(obs_para_list_nGenes == check_fct[["paralogs_list_nGenes"]])
    # stopifnot(obs_para_list_ratioTadsGenes == check_fct[["paralogs_list_ratioTadsGenes"]])
    
    # cat(paste0("stopifnot(obs_para_list_nTads == check_fct[[paralogs_list_nTads]])   \t=\t", 
    #            stopifnot(obs_para_list_nTads == check_fct[["paralogs_list_nTads"]]), "\n"))  
    # 
    # shuffle_chromoPartition <- function(domainDT, chrSize, preservePattern=TRUE, setSeed = F , seed = 42) {
    
    tad_dt <- read.delim(tad_file, header=FALSE, sep="\t", stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
    full_dt <- tad_dt
    tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
    tad_dt <- tad_dt[,c("chromo", "start", "end")]

    # iterate over number of permutations
    all_permut_results <- foreach(i_random=1:nRandom) %dopar% {
      
      
      # domain permutation for each chromosome
      rd_tad_dt <- foreach(chromo=c(unique(tad_dt$chromo)), .combine='rbind') %do% {
        chromo_tad_dt <- tad_dt[tad_dt$chromo == chromo,]
        chromo_rd_tad_dt <- shuffle_chromoPartition(domainDT=chromo_tad_dt, chrSize=max(chromo_tad_dt$end), preservePattern=FALSE, setSeed = FALSE ) 
        stopifnot(abs(sum(chromo_rd_tad_dt$end - chromo_rd_tad_dt$start) - sum(chromo_tad_dt$end - chromo_tad_dt$start)) < 1e-10)
        chromo_rd_tad_dt <- chromo_rd_tad_dt[order(chromo_rd_tad_dt$start, chromo_rd_tad_dt$end),]
        
        maxPos <- max(full_dt$end[full_dt$chromo == chromo])
        stopifnot(!is.na(maxPos) & is.numeric(maxPos))
        
        filled_chromo_rd_tad_dt <- fill_DT(dt=chromo_rd_tad_dt, chr_len=maxPos)
        
        stopifnot(filled_chromo_rd_tad_dt$region %in% c("TAD", "gap"))
        tmp_tad <- filled_chromo_rd_tad_dt[filled_chromo_rd_tad_dt$region == "TAD",]
        tmp_tad <- tmp_tad[order(tmp_tad$start, tmp_tad$end),]
        tmp_tad$region <- paste0(tmp_tad$chromo, "_TAD", 1:nrow(tmp_tad))
        
        tmp_gap <- filled_chromo_rd_tad_dt[filled_chromo_rd_tad_dt$region == "gap",]
        tmp_gap <- tmp_gap[order(tmp_gap$start, tmp_gap$end),]
        tmp_gap$region <- paste0(tmp_gap$chromo, "_BOUND", 1:nrow(tmp_gap))
        
        rd_chr_dt <- rbind(tmp_tad, tmp_gap)
        rd_chr_dt[order(rd_chr_dt$chromo, rd_chr_dt$start, rd_chr_dt$end),]
        
      }
      cat(paste0("... permut # ", i_random, "\t: end permut each chromo\n"))
      
      # do the gene2tad assignment based on TSS start
      # save(tad_dt, file="tad_dt.Rdata", version=2)
      # save(rd_tad_dt, file="rd_tad_dt.Rdata", version=2)
      # save(av_entrez, file="av_entrez.Rdata", version=2)
      # save(obs_g2t_dt, file="obs_g2t_dt.Rdata", version=2)
      
      # stopifnot(nrow(para_g2t_dt) > 0)
      av_gene=av_entrez[1]
      rd_g2t_dt <- foreach(av_gene = av_entrez, .combine='rbind') %do% {
        gchromo <- obs_g2t_dt$chromo[obs_g2t_dt$entrezID == av_gene]
        gstart <- obs_g2t_dt$start[obs_g2t_dt$entrezID == av_gene]
        gend <- obs_g2t_dt$end[obs_g2t_dt$entrezID == av_gene]
        stopifnot(!is.na(gstart) & is.numeric(gstart))
        stopifnot(!is.na(gend) & is.numeric(gend))
        
        # if(gstart > max(rd_tad_dt$end[rd_tad_dt$chromo == gchromo])) {
        #   rd_reg <- rd_tad_dt$region[rd_tad_dt$chromo == gchromo & rd_tad_dt$end == max(rd_tad_dt$end[rd_tad_dt$chromo == gchromo])]
        # } else {
        #   rd_reg <- rd_tad_dt$region[rd_tad_dt$chromo == gchromo & 
        #                                rd_tad_dt$start <= gstart & 
        #                                rd_tad_dt$end >= gstart]
        # }
        rd_reg <- rd_tad_dt$region[rd_tad_dt$chromo == gchromo & 
                                     rd_tad_dt$start <= gstart & 
                                     rd_tad_dt$end >= gstart]

        
        # cat(paste0("length(rd_reg)\t=", length(rd_reg), "\n"))
        
        stopifnot(length(rd_reg) == 1)
        
        data.frame(
          entrezID = as.character(av_gene),
          chromo = as.character(gchromo),
          region = as.character(rd_reg),
          stringsAsFactors = FALSE
        )
        
      }
      
      get_paralog_concord_with_tad(dt_para=entrez_paraDT,
                                                dt_g2t=rd_g2t_dt, 
                                                check_entrez=av_entrez)
    } # end-foreach iterating over permuts
    
    outFile <- file.path(outFolder, "all_permut_results.Rdata")
    save(all_permut_results, file=outFile, version=2)
    cat(paste0("... written: ", outFile, "\n"))
    
    list(obs_sameTAD_ratio=obs_sameTAD_ratio,
         obs_para_list_nTads=obs_para_list_nTads,
         obs_para_list_ratioTadsGenes =obs_para_list_ratioTadsGenes,
         permut_results = all_permut_results)
    
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
all_rd_ratio <- unlist(lapply(all_ds_results, function(sublist) lapply(sublist[["permut_results"]], function(x) x[["sameTAD_ratio"]])))

outFile <- file.path(outFolder, paste0("ratioParalogPairsInSameTAD.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth*1.2))  
plot_multiDens(
  list(
    obs_ratio = all_obs_ratio,
    permut_ratio = all_rd_ratio
  ),
  plotTit="ratio of paralog pairs in same TAD"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

all_obs_nTads <- unlist(lapply(all_ds_results, function(x)x[["obs_para_list_nTads"]]))
all_rd_nTads <- unlist(lapply(all_ds_results, function(sublist) lapply(sublist[["permut_results"]], function(x) x[["paralogs_list_nTads"]])))


outFile <- file.path(outFolder, paste0("ratioTadsGenes_paralogGroup.", plotType))
do.call(plotType, list(outFile, height = myHeight, width = myWidth*1.2))  
plot_multiDens(
  list(
    obs_nTads = all_obs_nTads,
    permut_nTads = all_rd_nTads
  ),
  plotTit="ratio of #TADs/#genes paralog groups"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


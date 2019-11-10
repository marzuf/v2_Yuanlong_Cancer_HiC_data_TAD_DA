startTime <- Sys.time()
cat(paste0("> Rscript paralogs_and_tads.R\n"))


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript paralogs_and_tads.R
buildTable <- TRUE


source("paralog_utils_fct.R")
nRandom <- 5
stopifnot(nRandom>1)
set.seed(10112019)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(pipFolder)
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

outFolder <- file.path("PARALOGS_AND_TADS")
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
  
  all_hicds = all_hicds[1]
  
  all_signif_genes_ic <- foreach(hicds = all_hicds) %dopar% {
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
          
              obs_para_list_nTads <- lengths(para_list_tads)
              hist(obs_para_list_nTads, breaks=1:max(obs_para_list_nTads))
    
    
    check_fct <- get_paralog_concord_with_tad(dt_para=entrez_paraDT,
                                             dt_g2t=obs_g2t_dt, 
                                             check_entrez=av_entrez)
      
    cat(paste0("check_fct[[sameTAD_ratio]] \t=\t", round(check_fct[["sameTAD_ratio"]], 4), "\n"))
    
    cat(paste0("stopifnot(obs_para_list_nTads == check_fct[[paralogs_list_nTads]])   \t=\t", 
               stopifnot(obs_para_list_nTads == check_fct[["paralogs_list_nTads"]]), "\n"))  
    
    # shuffle_chromoPartition <- function(domainDT, chrSize, preservePattern=TRUE, setSeed = F , seed = 42) {
    
    tad_dt <- read.delim(tad_file, header=FALSE, sep="\t", stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
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
        filled_chromo_rd_tad_dt <- fill_DT(dt=chromo_rd_tad_dt, chr_len=max(chromo_rd_tad_dt$end))
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
      # do the gene2tad assignment based on TSS start
      stopifnot(nrow(para_g2t_dt) > 0)
      av_gene=av_entrez[1]
      rd_g2t_dt <- foreach(av_gene = av_entrez, .combine='rbind') %do% {
        gchromo <- obs_g2t_dt$chromo[obs_g2t_dt$entrezID == av_gene]
        gstart <- obs_g2t_dt$start[obs_g2t_dt$entrezID == av_gene]
        gend <- obs_g2t_dt$end[obs_g2t_dt$entrezID == av_gene]
        stopifnot(!is.na(gstart) & is.numeric(gstart))
        stopifnot(!is.na(gend) & is.numeric(gend))
        
        rd_reg <- rd_tad_dt$region[rd_tad_dt$chromo == gchromo & 
                              rd_tad_dt$start <= gstart & 
                              rd_tad_dt$end >= gend]
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
    
  } # end-foreach hicds
} # END-IF BUILDTABLE
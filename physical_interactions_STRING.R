########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript physical_interactions_STRING.R\n"))

script_name <- "physical_interactions_STRING.R"


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript physical_interactions_STRING.R
buildTable <- TRUE

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

source("binding_interact_utils_fct.R")
nRandom <- 1000
# nRandom <- 1
stopifnot(nRandom>0)
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

outFolder <- file.path("PHYSICAL_INTERACTIONS_STRINGR")
dir.create(outFolder, recursive = TRUE)

########################################################################################################################################################################################


link_file <- "9606.protein.actions.v11.0.txt"
link_dt <- read.delim(link_file, header=TRUE, stringsAsFactors = FALSE)
link_dt$prot_id_a <- gsub("9606.", "", link_dt$item_id_a)
link_dt$prot_id_b <- gsub("9606.", "", link_dt$item_id_b)

id_file <- "ensembl_gene_prot.txt"
id_dt <- read.delim(id_file, sep="\t", stringsAsFactors=FALSE, header=TRUE)
id_dt <- na.omit(id_dt)
id_dt$EntrezGene.ID <- as.character(id_dt$EntrezGene.ID)

sum(link_dt$prot_id_a %in% id_dt$Protein.stable.ID)/nrow(link_dt) * 100
# 95%

sum(link_dt$prot_id_b %in% id_dt$Protein.stable.ID)/nrow(link_dt) * 100
# 95%

interact_mode <- "binding"

link_dt <- link_dt[link_dt$prot_id_a %in% id_dt$Protein.stable.ID & link_dt$prot_id_b %in% id_dt$Protein.stable.ID,]
nrow(link_dt)
link_dt <- link_dt[link_dt$mode == interact_mode,]

link_dt_entrezID <- left_join(link_dt[,c("prot_id_a", "prot_id_b")], id_dt[,c("Protein.stable.ID", "EntrezGene.ID")], by=c("prot_id_a"="Protein.stable.ID"))
colnames(link_dt_entrezID)[colnames(link_dt_entrezID) == "EntrezGene.ID"] <- "entrezID_a"
link_dt_entrezID <- left_join(link_dt_entrezID[,c("prot_id_a", "prot_id_b", "entrezID_a")], id_dt[,c("Protein.stable.ID", "EntrezGene.ID")], by=c("prot_id_b"="Protein.stable.ID"))
colnames(link_dt_entrezID)[colnames(link_dt_entrezID) == "EntrezGene.ID"] <- "entrezID_b"
head(link_dt_entrezID)



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
    
    tmp_link_dt_entrezID <- link_dt_entrezID[link_dt_entrezID$entrezID_a %in% obs_g2t_dt$entrezID & 
                                               link_dt_entrezID$entrezID_b %in% obs_g2t_dt$entrezID ,]
    
    av_entrez <- unique(c(tmp_link_dt_entrezID$entrezID_a, tmp_link_dt_entrezID$entrezID_b))
      
    
    x=get_binding_interact_with_tad(dt_para=tmp_link_dt_entrezID,
                                  dt_g2t=obs_g2t_dt, 
                                  check_entrez=av_entrez)
    
    
    
    
    interact_g2t_dt_tmp <- inner_join(tmp_link_dt_entrezID, obs_g2t_dt[,c("entrezID", "region")], by = c("entrezID_a" = "entrezID"))
    colnames(interact_g2t_dt_tmp)[colnames(interact_g2t_dt_tmp) == "region"] <- "region_a"
    interact_g2t_dt <- inner_join(interact_g2t_dt_tmp, obs_g2t_dt[,c("entrezID", "region")], by = c("entrezID_b" = "entrezID"))
    colnames(interact_g2t_dt)[colnames(interact_g2t_dt) == "region"] <- "region_b"
    nrow(interact_g2t_dt)
    # 158274
    stopifnot(is.character(interact_g2t_dt$entrezID_a))
    stopifnot(is.character(interact_g2t_dt$entrezID_b))
    stopifnot(setequal(av_entrez, unique(c(interact_g2t_dt$entrezID_a, interact_g2t_dt$entrezID_b))))
    
    obs_sameTAD_ratio <- sum(interact_g2t_dt$region_a == interact_g2t_dt$region_b)/nrow(interact_g2t_dt)

    stopifnot(x==obs_sameTAD_ratio)
    
    
    tad_dt <- read.delim(tad_file, header=FALSE, sep="\t", stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
    full_dt <- tad_dt
    tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
    tad_dt <- tad_dt[,c("chromo", "start", "end")]
    
    # iterate over number of permutations
    permut_sameTAD_ratios <- foreach(i_random=1:nRandom, .combine='c') %dopar% {
      
      
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
      
      
      av_gene=av_entrez[1]
      rd_g2t_dt <- foreach(av_gene = av_entrez, .combine='rbind') %do% {
        gchromo <- obs_g2t_dt$chromo[obs_g2t_dt$entrezID == av_gene]
        gstart <- obs_g2t_dt$start[obs_g2t_dt$entrezID == av_gene]
        gend <- obs_g2t_dt$end[obs_g2t_dt$entrezID == av_gene]
        stopifnot(!is.na(gstart) & is.numeric(gstart))
        stopifnot(!is.na(gend) & is.numeric(gend))
        
        
        rd_reg <- rd_tad_dt$region[rd_tad_dt$chromo == gchromo & 
                                     rd_tad_dt$start <= gstart & 
                                     rd_tad_dt$end >= gstart]
          
        stopifnot(length(rd_reg) == 1)
        
        data.frame(
          entrezID = as.character(av_gene),
          chromo = as.character(gchromo),
          region = as.character(rd_reg),
          stringsAsFactors = FALSE
        )
        
      }
      
      get_binding_interact_with_tad(dt_para=tmp_link_dt_entrezID,
                                   dt_g2t=rd_g2t_dt, 
                                   check_entrez=av_entrez)
    } # end-foreach iterating over permuts
    
    
    
    list(obs_sameTAD_ratio=obs_sameTAD_ratio,
         
         permut_sameTAD_ratios = permut_sameTAD_ratios)
    
  } # end-foreach hicds
  names(all_ds_results) <- all_hicds
  outFile <- file.path(outFolder, paste0(interact_mode, "_all_ds_results.Rdata"))
  save(all_ds_results, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
} else {# END-IF BUILDTABLE
  
  outFile <- file.path(outFolder, paste0(interact_mode, "_all_ds_results.Rdata"))
  all_ds_results <- get(load(outFile))
}

all_obs_ratio <- unlist(lapply(all_ds_results, function(x) x[["obs_sameTAD_ratio"]]))
all_rd_ratio <- unlist(lapply(all_ds_results, function(x) x[["permut_sameTAD_ratios"]]))

if(length(all_obs_ratio) == 1) {
  
  outFile <- file.path(outFolder, paste0("ratioInteractions_", interact_mode, "iInSameTAD.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth*1.2))  
  plot_multiDens(
    list(
      
      permut_ratio = all_rd_ratio
    ),
    plotTit="ratio of binding interactions in same TAD"
  )
  legend("topleft", legend=c("obs. ratio"), lty=2, col="red", bty="n")
  abline(v=all_obs_ratio, lty=2, col="red")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
} else{
  outFile <- file.path(outFolder, paste0("ratioInteractions_", interact_mode, "iInSameTAD.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth*1.2))  
  plot_multiDens(
    list(
      obs_ratio = all_obs_ratio,
      permut_ratio = all_rd_ratio
    ),
    plotTit="ratio of binding interactions in same TAD"
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}





##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


conserved_dt <- read.delim("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2/conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3_df.txt",
                           stringsAsFactors = FALSE, header=T)

conserved_dt$conserved_exprds <- sapply(conserved_dt$corresp_tads, function(x) paste0(sort(unique(sapply(strsplit(x = x, split=","), function(k) basename(dirname(k))))), collapse=","))


conserved_dt$nConserved_dt <- sapply(conserved_dt$corresp_tads, function(x) length(unique(sapply(strsplit(x = x, split=","), function(k) k))))

conserved_dt$nConserved_exprds <- sapply(conserved_dt$corresp_tads, function(x) length(unique(sapply(strsplit(x = x, split=","), function(k) basename(dirname(k))))))

 View(conserved_dt[order(conserved_dt$nConserved_exprds, decreasing = TRUE),])
 
 View(conserved_dt[order(conserved_dt$nConserved_dt, decreasing = TRUE),])






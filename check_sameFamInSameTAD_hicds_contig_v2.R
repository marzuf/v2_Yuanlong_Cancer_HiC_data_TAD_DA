
startTime <- Sys.time()

### !! v2 version: 
# - contiguty i.e. gene rank based on all genes (not only genes for which I have annotation)
# - sample genes among all genes, not only those for which I have annotation

# Rscript check_sameFamInSameTAD_hicds_contig_v2.R # RUN POSITRON

# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity

# no need to have exprds

set.seed(20200903) # 

require(doMC)
require(foreach)
registerDoMC(40)
require(reshape2)
require(igraph)

runFolder <- "."
pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")
familyVar <- "hgnc_family_short"

corMethod <- "pearson"
familyData <- "hgnc_family_short"

nRandom <- 100
# nRandom=5

logOffset <- 0.01


plotType <- "svg"
myHeight <- 5
myWidth <- 7

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

gff_dt$true_start <- ifelse(as.character(gff_dt$strand) == "+", gff_dt$start, gff_dt$end)
stopifnot(is.numeric(gff_dt$true_start))
gff_dt <- gff_dt[order(gff_dt$chromo, gff_dt$true_start, gff_dt$end),]



all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
# all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))


# all_ds <- unlist(sapply(all_hicds, function(x) file.path(x, list.files(file.path(pipFolder, x)))), use.names = FALSE)

# hicds = "Barutcu_MCF-10A_40kb"
# exprds = "TCGAbrca_lum_bas"
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]

buildData <- TRUE

aggFunc <- "mean"

### FOR CHECK !!!
setDir <- "/media/electron"
setDir <- ""
hgnc_geneFamilyFile <- file.path(setDir, "/mnt/ed4/marie/family_data_2/hgnc_entrez_family.txt")
hgnc_geneFamilyDT <- read.delim(hgnc_geneFamilyFile, col.names=c("entrezID", "family"), header = F, stringsAsFactors = F)
hgnc_geneFamilyDT$entrezID <- as.character(hgnc_geneFamilyDT$entrezID)
hgnc_geneFamilyDT$family_short <- unlist(sapply(hgnc_geneFamilyDT$family, function(x) strsplit(x, "\\|")[[1]][1] ))
# any(duplicated(hgnc_geneFamilyDT$entrezID))


# ds=all_ds[1]

# all_hicds=all_hicds[1]

if(buildData) {
  
  
  all_ds_results <- foreach(hicds = all_hicds) %do%{
    
    
    cat(paste0("... start: ", hicds,"\n"))
    
    
    #       ### => CHANGED FOR THE TISSUE DATA TO USE TISSUE SPECIFIC FAMILY FILES !!!
    #       # inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
          inFoldFamily <- file.path("PREP_GENE_FAMILIES_TAD_DATA", hicds)
          familyData2 <- "hgnc"
          familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData2, "_entrezID_family_TAD_DT.Rdata")))))
          familyDT$entrezID <- as.character(familyDT$entrezID)
    
    g2t_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt")
    g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    tad_g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    stopifnot(!duplicated(g2t_dt$entrezID))
    entrezIDchromo <- setNames(g2t_dt$chromo, g2t_dt$entrezID)
    
    entrezIDregion <- setNames(tad_g2t_dt$chromo, tad_g2t_dt$entrezID)
    
    tad_entrezID <- tad_g2t_dt$entrezID  ### v2 => I will sample from them !!!
    stopifnot(tad_entrezID %in% names(entrezIDchromo))
    tad_entrezIDchromo <- entrezIDchromo[names(entrezIDchromo) %in% tad_entrezID] # used for the sampling
    
    
    sameTADfile <- file.path("CREATE_SAME_TAD_SORTNODUP", hicds, "all_TAD_pairs.Rdata")
    stopifnot(file.exists(sameTADfile))
    sameTAD_dt <- get(load(sameTADfile))
    
    #stopifnot(tad_entrezID %in% sameTAD_dt$gene1 | tad_entrezID %in% sameTAD_dt$gene2 )
    # => will not be true for genes from TAD that contain only 1 gene
    big_tads <- names(table(tad_g2t_dt$region)) [as.numeric(table(tad_g2t_dt$region)) > 1]
    big_tad_g2t_dt <- tad_g2t_dt[tad_g2t_dt$region %in% big_tads,]
    stopifnot(tad_entrezID[tad_entrezID %in% big_tad_g2t_dt$entrezID] %in% sameTAD_dt$gene1 |
                tad_entrezID[tad_entrezID %in% big_tad_g2t_dt$entrezID] %in% sameTAD_dt$gene2 )
    
    
    # ADDED 08.01.19 to accommodate updated family file
    sameFamFolder <- file.path("CREATE_SAME_FAMILY_SORTNODUP", hicds)
    # checking the file comes after (iterating over family and family_short)
    stopifnot(dir.exists(sameFamFolder))
    sameFamFile <- file.path(sameFamFolder, paste0(familyData, "_all_family_pairs.Rdata")) # at least this one should exist !
    stopifnot(file.exists(sameFamFile))
    sameFam_dt <- get(load(sameFamFile))
    
    stopifnot(sameTAD_dt$gene1 <= sameTAD_dt$gene2)
    stopifnot(sameFam_dt$gene1 <= sameFam_dt$gene2)
    stopifnot(sameFam_dt$family %in% familyDT[,familyData])
    
    stopifnot(sameFam_dt$gene1 %in% g2t_dt$entrezID) 
    stopifnot(sameFam_dt$gene2 %in% g2t_dt$entrezID)
    stopifnot(sameFam_dt$gene1 %in% tad_g2t_dt$entrezID)
    stopifnot(sameFam_dt$gene2 %in% tad_g2t_dt$entrezID)
    
    
    sameFamSameTAD_dt <- merge(sameFam_dt, sameTAD_dt, all.x=TRUE, all.y=FALSE, by=c("gene1", "gene2"))
    
    all_genes <- familyDT$entrezID
    stopifnot(all_genes %in% names(entrezIDchromo))
    family_entrezIDchromo <- entrezIDchromo[names(entrezIDchromo) %in% all_genes] # used for the sampling - not in V2 !!!
    stopifnot(names(family_entrezIDchromo) %in% tad_g2t_dt$entrezID)
    rm("family_entrezIDchromo") # not used in v2 !!!
    
    all_fams <- unique(sameFamSameTAD_dt$family)
    # all_fams=all_fams[1]
    
    ## !!! added here for gene rank
    # hicds_gff_dt <- gff_dt[gff_dt$entrezID %in% names(family_entrezIDchromo),]
    # ADAPTED 19.10.20 V2: CONTIGUITY BASED ON ALL GENES ASSIGNED TO A TAD
    stopifnot(tad_entrezID %in% gff_dt$entrezID)
    hicds_gff_dt <- gff_dt[gff_dt$entrezID %in% tad_entrezID,]
    hicds_gff_dt <- hicds_gff_dt[order(hicds_gff_dt$chromo, hicds_gff_dt$true_start, hicds_gff_dt$end),]

      
    # to avoid iterating over chromo -> when switching chromo, add a step of 2
    hicds_gff_dt$chromo <- as.character(hicds_gff_dt$chromo)
    
    # 17.10.2020: changed -> not needed, because i do this chromosome by chromosome
    # otherwise i could not select adjacent positions in the ranks
    # hicds_gff_dt$gene_rank_offset <- sapply(hicds_gff_dt$chromo, function(x) which(unique(hicds_gff_dt$chromo) == x))    
    # hicds_gff_dt$gene_rank_offset <- hicds_gff_dt$gene_rank_offset - 1
    # hicds_gff_dt$gene_rank_init <- 1:nrow(hicds_gff_dt)
    # hicds_gff_dt$gene_rank <- hicds_gff_dt$gene_rank_offset + hicds_gff_dt$gene_rank_init
    # finally i dont need the break for the chromo because I sample by chromo
    # hicds_gff_dt[9121:9123,]
    # hicds_gff_dt[12078:12080,]
    # stopifnot(sum(diff(hicds_gff_dt$gene_rank) != 1) == length(unique(hicds_gff_dt$chromo))-1)
    # stopifnot(sum(diff(hicds_gff_dt$gene_rank) == 2) == length(unique(hicds_gff_dt$chromo))-1)
    hicds_gff_dt$gene_rank <- 1:nrow(hicds_gff_dt)
    
    entrezID_geneRanks <- setNames(hicds_gff_dt$gene_rank, hicds_gff_dt$entrezID)
    
    # stopifnot(setequal(names(family_entrezIDchromo), names(entrezID_geneRanks)))
    # stopifnot(names(family_entrezIDchromo) %in% hicds_gff_dt$entrezID)
    # v2 19.10.2020
    stopifnot(setequal(names(tad_entrezIDchromo), names(entrezID_geneRanks)))
    stopifnot(names(tad_entrezIDchromo) %in% hicds_gff_dt$entrezID)
    
    
    ##>> iterate here over families
    i_fam=1
    all_fam_results <- foreach(i_fam = 1:length(all_fams)) %dopar% {
      # all_fam_results <- foreach(i_fam = 1) %dopar% {
      
      fam <- all_fams[i_fam]
      
      cat(paste0("... ", hicds, " - start fam: ", i_fam, "/", length(all_fams), "\n"))
      
      fam_dt <- sameFamSameTAD_dt[sameFamSameTAD_dt$family == fam,]
      
      sameTAD_fam_dt <- fam_dt[!is.na(fam_dt$region),]
      
      nbrSingletons <- length(setdiff(c(fam_dt$gene1,fam_dt$gene2), c(sameTAD_fam_dt$gene1, sameTAD_fam_dt$gene2)))
      
      nbrGenes <- length(unique(c(fam_dt$gene1, fam_dt$gene2)))
      
      stopifnot(nbrSingletons + length(unique(c(sameTAD_fam_dt$gene1, sameTAD_fam_dt$gene2))) == nbrGenes)
      
      stopifnot(nbrSingletons <= nbrGenes)
      
      fam_net <- graph_from_data_frame(d=sameTAD_fam_dt[,c("gene1", "gene2")], directed=F) 
      nbrEdges <- gsize(fam_net)
      stopifnot(nbrEdges == nrow(sameTAD_fam_dt))
      
      nbrComponents <- components(fam_net)$no
      stopifnot(nbrComponents == length(unique(sameTAD_fam_dt$region)))
      
      true_genes <- unique(c(fam_dt$gene1, fam_dt$gene2))
      stopifnot(true_genes %in% names(entrezIDchromo))
      true_chromos <- entrezIDchromo[true_genes]
      true_countChromos <- table(true_chromos)
      
      # script0_name <- "0_prepGeneData"
      # pipeline_geneList <- get(load(file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")))
      # all(familyDT$entrezID %in% pipeline_geneList) # FALSE
      # all(sameTAD_dt$gene1 %in% pipeline_geneList) # FALSE
      
      # FOR THE MOMENT -> DONT DO THIS AT CHROMOSOME LEVEL
      ### added 17.10.20 -> retrieve how many contiguous genes
      stopifnot(true_genes %in% names(entrezID_geneRanks))
      
      ## should be done chromo by chromo !!!!
      # true_ranks <- sort(entrezID_geneRanks[true_genes])
      # stopifnot(length(true_ranks) == length(true_genes))
      # names(true_ranks) <- NULL
      # cont_gene_ranks <- diff(true_ranks)
      # stopifnot(cont_gene_ranks > 0)
      # rle_cont_genes <- rle(cont_gene_ranks)
      # contig_tosample <- rle_cont_genes$lengths[rle_cont_genes$values==1]+1
      # nbr_notcontig_tosample <- length(true_genes) - sum(contig_tosample)
      
      
      
      # save(sameFamSameTAD_dt, file="sameFamSameTAD_dt.Rdata", version=2)
      # load( file="sameFamSameTAD_dt.Rdata")
      # save(sameFamSameTAD_dt, file="sameFamSameTAD_dt.Rdata", version=2)
      # load( file="sameFamSameTAD_dt.Rdata")
      # save(sameTAD_fam_dt, file="sameTAD_fam_dt.Rdata", version=2)
      # load(file="sameTAD_fam_dt.Rdata")
      # save(family_entrezIDchromo, file="family_entrezIDchromo.Rdata", version=2)
      # load( file="family_entrezIDchromo.Rdata")
      # save(fam_dt, file="fam_dt.Rdata", version=2)
      # load(file="fam_dt.Rdata")
      # stop("-ok\n")
      
      random_results_dt <- foreach(i = 1:nRandom, .combine='cbind') %do% {
        
        # # V0: SAMPLE CHROMO BY CHROMO
        # sample_genes <- foreach(chromo = names(true_countChromos), .combine='c') %do% {
        #   sample(names(family_entrezIDchromo)[family_entrezIDchromo == chromo], 
        #          size=true_countChromos[chromo], 
        #          replace=FALSE)
        # }
        # stopifnot(!duplicated(sample_genes))
        # V2: SAMPLE CONTIGUOUS - FOR THE MOMENT, NOT CHROMOSOME-WISE
        # !!! changed: i need to it chromosome wise because start_pos:start_pos+x could go over chromosomes otherwise
        
        sample_genes <- foreach(chromo = names(true_countChromos), .combine='c') %do% {
          # v2: change here 19.10.2020 -> sample from all genes
          # all_chromo_genes <- names(family_entrezIDchromo)[family_entrezIDchromo == chromo]
          all_chromo_genes <- names(tad_entrezIDchromo)[tad_entrezIDchromo == chromo]
           
           true_chromo_genes <- true_genes[true_genes %in% all_chromo_genes]
           stopifnot(length(true_chromo_genes) > 0)
           stopifnot(all_chromo_genes %in% names(entrezID_geneRanks))
           stopifnot(true_chromo_genes %in% names(entrezID_geneRanks))
           
           # bag in which I sample:
           chromo_entrezID_geneRanks <- entrezID_geneRanks[names(entrezID_geneRanks) %in% all_chromo_genes]
           
           
           
           true_ranks <- sort(chromo_entrezID_geneRanks[true_chromo_genes])
           stopifnot(length(true_ranks) == length(true_chromo_genes))
           names(true_ranks) <- NULL
           cont_gene_ranks <- diff(true_ranks)
           stopifnot(cont_gene_ranks > 0)
           rle_cont_genes <- rle(cont_gene_ranks)
           contig_tosample <- rle_cont_genes$lengths[rle_cont_genes$values==1]+1  # this works even if empty : if c() + 1 = numeric(0)
           nbr_notcontig_tosample <- length(true_chromo_genes) - sum(contig_tosample)
           
           
           # contig_sampled <- unlist( sapply(contig_tosample, function(x) {
           #   # sample the start index
           #   start_pos <- sample(x = 1:(length(chromo_entrezID_geneRanks)-x), size = 1)
           #   samp_genes <- names(chromo_entrezID_geneRanks)[start_pos:(start_pos+x-1)]
           #   stopifnot(!is.na(samp_genes))
           #   samp_genes
           # }))
           # I CANNOT ENSURE IT BECAUSE THE WAY I CONSTRUCT
           # CONTIG SAMPLED I DO NOT ENSURE
           # stopifnot(!duplicated(sample_genes))
           # => need a for-loop because I should ensure not duplicated !!!
           if(length(contig_tosample) == 0) {
             contig_sampled <- c()
           } else {
             tmp_entrezID_geneRanks <- chromo_entrezID_geneRanks
             contig_sampled <- list()
             for(k in 1:length(contig_tosample)) {
               nToSamp <- contig_tosample[k]
               stopifnot(length(tmp_entrezID_geneRanks) >= nToSamp)
               # for a vector of 3, if nSamp=2, maxIdx=2 (3-2+1)
               # for a vector of 5, if nSamp=3, maxIdx=3 (5-3+1)
               start_pos <- sample(x = 1:(length(tmp_entrezID_geneRanks)-nToSamp+1), size = 1)
               # if startIdx=2, if nSamp=2, sample 2:(2+2-1)
               samp_genes <- names(tmp_entrezID_geneRanks)[start_pos:(start_pos+nToSamp-1)]
               stopifnot(!is.na(samp_genes))
               tmp_entrezID_geneRanks <- tmp_entrezID_geneRanks[! names(tmp_entrezID_geneRanks) %in% samp_genes]
               contig_sampled[[k]] <- samp_genes
             }
             contig_sampled <- unlist(contig_sampled)
             stopifnot(!duplicated(contig_sampled))
             stopifnot(length(tmp_entrezID_geneRanks) == length(chromo_entrezID_geneRanks) - sum(contig_tosample))
             stopifnot(length(contig_sampled) == sum(contig_tosample))
             stopifnot(contig_sampled %in% names(chromo_entrezID_geneRanks))
             
           }
           if(nbr_notcontig_tosample > 0) {
             stopifnot(length(contig_sampled) == sum(contig_tosample))
             notcontig_sampled <- sample(
               x = names(chromo_entrezID_geneRanks)[!names(chromo_entrezID_geneRanks) %in% contig_sampled], # DISCARD THOSE ALREADY SAMPLED
               size = nbr_notcontig_tosample
             )
             stopifnot(length(notcontig_sampled) == nbr_notcontig_tosample)
             stopifnot(!duplicated(notcontig_sampled))
             
           } else {
             notcontig_sampled <- c()
           }
           chromo_sample_genes <- c(contig_sampled, notcontig_sampled)
           stopifnot(!duplicated(chromo_sample_genes))
           stopifnot(chromo_sample_genes %in% names(chromo_entrezID_geneRanks))
           stopifnot(chromo_sample_genes %in% all_chromo_genes)
           chromo_sample_genes
        }
        
        stopifnot(length(sample_genes) == length(true_genes))
        
        sample_dt <- as.data.frame(t(combn(x=sample_genes, m = 2)))
        colnames(sample_dt) <- c("gene1_tmp", "gene2_tmp")
        
        sample_dt$gene1_tmp <- as.character(sample_dt$gene1_tmp)
        sample_dt$gene2_tmp <- as.character(sample_dt$gene2_tmp)
        sample_dt$gene1 <-as.character(pmin(sample_dt$gene1_tmp, sample_dt$gene2_tmp))
        sample_dt$gene2 <-as.character(pmax(sample_dt$gene1_tmp, sample_dt$gene2_tmp))
        stopifnot(sample_dt$gene1 == sample_dt$gene1_tmp | sample_dt$gene1 == sample_dt$gene2_tmp)
        stopifnot(sample_dt$gene2 == sample_dt$gene1_tmp | sample_dt$gene2 == sample_dt$gene2_tmp)
        sample_dt$gene1_tmp <- sample_dt$gene2_tmp <- NULL
        stopifnot(sample_dt$gene1 < sample_dt$gene2)
        
        random_sameTAD_fam_dt <- merge(sample_dt, sameTAD_dt, all.x=T, all.y=F, by=c("gene1", "gene2"))
        random_sameTAD_fam_dt <- random_sameTAD_fam_dt[!is.na(random_sameTAD_fam_dt$region),]
        
        random_sameTAD_fam_dt$gene1 <- as.character(random_sameTAD_fam_dt$gene1)
        stopifnot(random_sameTAD_fam_dt$gene1 %in% names(entrezIDregion))
        random_sameTAD_fam_dt$gene2 <- as.character(random_sameTAD_fam_dt$gene2)
        stopifnot(random_sameTAD_fam_dt$gene2 %in% names(entrezIDregion))
        
        all_regions_gene1 <- entrezIDregion[random_sameTAD_fam_dt$gene1]
        all_regions_gene2 <- entrezIDregion[random_sameTAD_fam_dt$gene2]
        
        stopifnot(all_regions_gene1 == all_regions_gene2)
        
        return(NULL)
      } # end permut
      
      return(NULL)
    } # end families

    return(NULL)
  } # end datasets
}

cat("# check ok\n")
cat(paste0("*** DONE\n", startTime, "\n", Sys.time(), "\n"))


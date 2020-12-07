
startTime <- Sys.time()

# Rscript sameFamInSameTAD_hicds_contig_disc_v2.R  # RUN ELECTRON

### !! v2 version: 
# - contiguty i.e. gene rank based on all genes (not only genes for which I have annotation)
# - sample genes among all genes, not only those for which I have annotation
# - discard if 0.8 of the genes belong to same family  

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

maxSameFamRatioThresh <- 0.8
# keep sample if: all(fams_ratio < maxSameFamRatioThresh )

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


outFolder <- file.path("SAMEFAMINSAMETAD_HICDS_CONTIG_DISC_V2")
dir.create(outFolder, recursive = TRUE)

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

buildData <- FALSE

aggFunc <- "mean"
logOffset <- 0.01

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
    fam2entrezID <- setNames(familyDT[,paste0(familyData)], familyDT$entrezID)
    
    g2t_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt")
    g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    tad_g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    stopifnot(!duplicated(g2t_dt$entrezID))
    entrezIDchromo <- setNames(g2t_dt$chromo, g2t_dt$entrezID)
    
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
    # all_fams=all_fams[1:5]
    
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
          cont_gene_ranks <- diff(true_ranks) # remark 7.12: I think I could have directly do rle(...$familyID) without using diff()
          stopifnot(cont_gene_ranks > 0)
          rle_cont_genes <- rle(cont_gene_ranks)
          contig_tosample <- rle_cont_genes$lengths[rle_cont_genes$values==1]+1
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
              check_inf <- 0
              while(TRUE) {
                start_pos <- sample(x = 1:(length(tmp_entrezID_geneRanks)-nToSamp+1), size = 1)
                # if startIdx=2, if nSamp=2, sample 2:(2+2-1)
                samp_genes <- names(tmp_entrezID_geneRanks)[start_pos:(start_pos+nToSamp-1)]
                stopifnot(!is.na(samp_genes))
                ### for _disc -> discard if same family added 18.10.2020
                # stopifnot(samp_genes %in% names(fam2entrezID)) # v2 not true since sample from all genes
                # v2: 19.10.20
                stopifnot(samp_genes %in% names(tad_entrezIDchromo))
                stopifnot(samp_genes %in% tad_entrezID)
                # samp_fams <- as.character(fam2entrezID[samp_genes])  # v2 not possible since sample from all genes
                samp_fams <- fam2entrezID[names(fam2entrezID) %in% samp_genes]
                save(samp_genes, file="samp_genes.Rdata", version=2)
                save(samp_fams, file="samp_fams.Rdata", version=2)
                # v2: no family with > 80% of the genes
                fams_ratio <- as.numeric(table(samp_fams))/length(samp_genes) # and not length samp_fams since sample from ALL genes
                # NB: it never happens that length(samp_genes) == 1
                # otherwise not contiguous ;) ; checked in LOOK_SAMEFAM_CONTIG - stopifnot added 19.10
                stopifnot(length(samp_genes) > 1)
                # if(length(samp_genes) == 1 | length(unique(samp_fams)) > 1) {
                # changed for v2
                if(length(samp_genes) == 1 | all(fams_ratio < maxSameFamRatioThresh )) {  
                  tmp_entrezID_geneRanks <- tmp_entrezID_geneRanks[! names(tmp_entrezID_geneRanks) %in% samp_genes]
                  contig_sampled[[k]] <- samp_genes
                  break
                }
                check_inf <- check_inf+1
                stopifnot(check_inf < 100)
              }
            }
            stopifnot(!unlist(lapply(contig_sampled, is.null)))
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
        
        # NB 19.10 => would be faster to do the combn inside chromo because they will never be in same TAD if from different chromo...
        
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
        
        random_sameTAD_fam_dt <- merge(sample_dt, sameTAD_dt, all.x=T, all.y=F, by=c("gene1", "gene2"))
        random_sameTAD_fam_dt <- random_sameTAD_fam_dt[!is.na(random_sameTAD_fam_dt$region),]
        
        random_nbrSingletons <- nbrGenes -length(unique(c(random_sameTAD_fam_dt$gene1, random_sameTAD_fam_dt$gene2)))
        
        random_fam_net <- graph_from_data_frame(d=random_sameTAD_fam_dt[,c("gene1", "gene2")], directed=F) 
        random_nbrEdges <- gsize(random_fam_net)
        stopifnot(random_nbrEdges == nrow(random_sameTAD_fam_dt))
        
        random_nbrComponents <- components(random_fam_net)$no
        stopifnot(random_nbrComponents == length(unique(random_sameTAD_fam_dt$region)))
        
        c(random_nbrSingletons=random_nbrSingletons, random_nbrEdges=random_nbrEdges, random_nbrComponents=random_nbrComponents)
      } # end permut
      
      list(
        nbrGenes=nbrGenes,
        nbrSingletons=nbrSingletons,
        nbrEdges=nbrEdges,
        nbrComponents=nbrComponents,
        random_results_dt=random_results_dt
      )
    } # end families
    names(all_fam_results) <- all_fams
    outFile <- file.path(outFolder, paste0(hicds, "_", "all_fam_results.Rdata"))
    save(all_fam_results, file=outFile, version=2)
    all_fam_results
  } # end datasets
  names(all_ds_results) <- all_hicds
  
  outFile <- file.path(outFolder, "all_ds_results.Rdata")
  save(all_ds_results, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_ds_results.Rdata")
  all_ds_results <- get(load(outFile))
}

all_obs_nbrEdges <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["nbrEdges"]])))
all_obs_nbrCpts <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["nbrComponents"]])))
all_obs_nbrSingletons <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["nbrSingletons"]])))

length(all_obs_nbrEdges)
sum(all_obs_nbrEdges == 0)

length(all_obs_nbrCpts)
sum(all_obs_nbrCpts == 0)


all_random_nbrEdges <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["random_results_dt"]]['random_nbrEdges',])))
all_random_meanNbrEdges <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) mean(x[["random_results_dt"]]['random_nbrEdges',]))))
stopifnot(length(all_obs_nbrEdges) == length(all_random_meanNbrEdges))

length(all_random_nbrEdges)
sum(all_random_nbrEdges == 0)


all_random_nbrCpts <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["random_results_dt"]]['random_nbrComponents',])))
all_random_meanNbrCpts <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) mean(x[["random_results_dt"]]['random_nbrComponents',]))))
stopifnot(length(all_obs_nbrCpts) == length(all_random_meanNbrCpts))

all_random_nbrSingletons <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["random_results_dt"]]['random_nbrSingletons',])))
all_random_meanNbrSingletons <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) mean(x[["random_results_dt"]]['random_nbrSingletons',]))))
stopifnot(length(all_obs_nbrSingletons) == length(all_random_meanNbrSingletons))

nPermut <- unique(unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) ncol(x[["random_results_dt"]])))))
stopifnot(length(nPermut) == 1)


length(all_random_nbrCpts)
sum(all_random_nbrCpts == 0)

source("../Cancer_HiC_data_TAD_DA/utils_plot_fcts.R")

### ITERATE TO DO THE SAME FOR EACH OF THE VARIABLE!!!

plotCex <- 1.2

all_vars <- c("Edges", "Cpts", "Singletons")
all_vars <- c("Edges")
curr_var=all_vars[1]
for(curr_var in all_vars) {
  
  
  # aggregate by family
  nbr_obs_dt <- data.frame(family_lab=names(eval(parse(text=paste0("all_obs_nbr", curr_var)))), 
                           nbr=as.numeric(eval(parse(text=paste0("all_obs_nbr", curr_var)))), stringsAsFactors = FALSE)
  nbr_obs_dt$family <- gsub(".+?\\.(.+)", "\\1", nbr_obs_dt$family_lab)  ### CORRECT HERE THERE DOTS IN THE FAMILY NAMES SHOULD USE INTERROGATION MARK 25.09.2020
  ### !!!! ADDED 25.09 IMPORTANT CHECK
  stopifnot(nbr_obs_dt$family %in% hgnc_geneFamilyDT$family_short)
  
  stopifnot(nbr_obs_dt$family != "")
  stopifnot(!is.na(nbr_obs_dt$family))
  
  nbr_rd_dt <- data.frame(family_lab=names(eval(parse(text=paste0("all_random_nbr", curr_var)))), 
                          nbr=as.numeric(eval(parse(text=paste0("all_random_nbr", curr_var)))), 
                          stringsAsFactors = FALSE)
  nbr_rd_dt$family <- gsub(".+?\\.(.+)\\.result.+", "\\1", nbr_rd_dt$family_lab)  ### CORRECT HERE THERE DOTS IN THE FAMILY NAMES SHOULD USE INTERROGATION MARK 25.09.2020
  stopifnot(nbr_rd_dt$family != "")
  stopifnot(!is.na(nbr_rd_dt$family))
  ### !!!! ADDED 25.09 IMPORTANT CHECK
  stopifnot(nbr_rd_dt$family %in% hgnc_geneFamilyDT$family_short)
  
  stopifnot(setequal(nbr_obs_dt$family, nbr_rd_dt$family))
  
  
  agg_obs_dt <- aggregate(nbr~family, data=nbr_obs_dt, FUN=aggFunc)
  agg_rd_dt <- aggregate(nbr~family, data=nbr_rd_dt, FUN=aggFunc)
  
  
  plot_dt <- merge(agg_obs_dt, agg_rd_dt, by="family", all=T, suffixes=c("_obs", "_rd"))
  stopifnot(!is.na(plot_dt))
  
  plotTit <- paste0("# of ", curr_var, " by family")
  subTit <- paste0("aggreg. func = ", aggFunc, "; # families = ", length(unique(plot_dt$family)))
  
  plot_dt$nbr_obs_log10 <- log10(plot_dt$nbr_obs + logOffset)
  plot_dt$nbr_rd_log10 <- log10(plot_dt$nbr_rd + logOffset)
  
  save(plot_dt, file=file.path(outFolder, "sameFam_nbrEdges_contig_disc_plot_dt_scatterplot.Rdata"), version=2)
  cat(paste0("... written: ", file.path(outFolder, "sameFam_nbrEdges_contig_disc_plot_dt_scatterplot.Rdata"),  "\n"))
  
  
  my_x <- plot_dt$nbr_obs_log10
  my_y <- plot_dt$nbr_rd_log10
  
  outFile <- file.path(outFolder,  paste0("allDS_nbr", curr_var, "_byFam_scatterplot_", gsub("\\.", "", logOffset), ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myHeight))
  plot(x=my_x,y=my_y, main=plotTit, 
       pch=16, cex=0.7,
       cex.axis=plotCex,
       cex.lab=plotCex,
       cex.main=plotCex,
       xlab=paste0("# ", curr_var, " obs. [log10(+", logOffset, ")]"), 
       ylab=paste0("# ", curr_var, " rd. [log10(+", logOffset, ")]"))
  curve(1*x, col="grey", add=TRUE)
  addCorr(x=my_x, y=my_y, legPos="topleft", bty="n")
  mtext(side=3, text = subTit, font=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile,  "\n"))
  
  my_x <- plot_dt$nbr_obs
  my_y <- plot_dt$nbr_rd
  
  outFile <- file.path(outFolder,  paste0("allDS_nbr", curr_var, "_byFam_scatterplot_notLog.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myHeight))
  plot(x=my_x,y=my_y, main=plotTit, 
       pch=16, cex=0.7,
       cex.axis=plotCex,
       cex.lab=plotCex,
       cex.main=plotCex,
       xlab=paste0("# ", curr_var, " obs."), 
       ylab=paste0("# ", curr_var, " rd."))
  curve(1*x, col="grey", add=TRUE)
  addCorr(x=my_x, y=my_y, legPos="topleft", bty="n")
  mtext(side=3, text = subTit, font=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile,  "\n"))
  
}

cat(paste0("*** DONE\n", startTime, "\n", Sys.time(), "\n"))


# # all_vars <- c("Edges", "Cpts", "Singletons")
# all_vars <- c("Edges")
# curr_var=all_vars[1]
# for(curr_var in all_vars) {
#   
#   plotTit <- paste0("nbr", curr_var)
#   subTit <- paste0("all DS - n =", length(all_ds_results))
#   
#   
#   plot_list <- list(log10(logOffset+eval(parse(text = paste0("all_obs_nbr", curr_var)))),
#                     log10(logOffset+eval(parse(text = paste0("all_random_nbr", curr_var)))))
#   
#   
#   
#   names(plot_list) <- c(paste0("all_obs_nbr",curr_var), paste0("all_random_nbr", curr_var) )
#   
#   outFile <- file.path(outFolder,  paste0("allDS_nbr", curr_var, "_density_", gsub("\\.", "", logOffset), ".", plotType))
#   do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#   plot_multiDens(
#     plot_list,
#     plotTit = plotTit, 
#     my_xlab = paste0("# ", curr_var, " [log10(+", logOffset, ")]")
#   )
#   mtext(side=3, text = subTit, font=3)
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile,  "\n"))
#   
#   
#   
#   plotTit <- paste0("nbr", curr_var, "[log10(+", logOffset, ")]")
#   subTit <- paste0("all DS - n =", length(all_ds_results))
#   
#   xlab <- "observed"
#   ylab <- paste0("mean permut. (n=", nPermut, ")" )
#   
#   # outFile <- file.path(outFolder,  paste0("allDS_nbr", curr_var, "_meanPermut_vs_obs_densplot.", "svg"))
#   outFile <- file.path(outFolder,  paste0("allDS_nbr", curr_var, "_meanPermut_vs_obs_densplot_", gsub("\\.", "", logOffset), ".", "svg"))
#   do.call("svg", list(outFile, height=7, width=7))
#   
#   
#   densplot(x=log10(eval(parse(text=paste0("all_obs_nbr", curr_var))) +logOffset),
#            y=log10(eval(parse(text=paste0("all_random_meanNbr", curr_var))) +logOffset),
#            xlab= xlab,ylab=ylab, main=plotTit
#   )
#   mtext(side=3, text = subTit, font=3)
#   curve(1*x, lty=1, col="darkgrey", add = T)
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile,  "\n"))
#   
#   
# }
# 

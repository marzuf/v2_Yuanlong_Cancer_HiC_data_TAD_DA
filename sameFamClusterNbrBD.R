
### !! v2 version: 
# - contiguty i.e. gene rank based on all genes (not only genes for which I have annotation)
# - sample genes among all genes, not only those for which I have annotation

# Rscript sameFamClusterNbrBD.R # RUN POSITRON

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

all_entrezIDstarts <- setNames(gff_dt$true_start, gff_dt$entrezID) #


outFolder <- file.path("SAMEFAMCLUSTERNBRBD")
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

# all_hicds=all_hicds[2]

if(buildData) {
  
  
  all_ds_results <- foreach(hicds = all_hicds) %do%{
    
    
    cat(paste0("... start: ", hicds,"\n"))
    
    
    #       ### => CHANGED FOR THE TISSUE DATA TO USE TISSUE SPECIFIC FAMILY FILES !!!
    #       # inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
    inFoldFamily <- file.path("PREP_GENE_FAMILIES_TAD_DATA", hicds)
    familyData2 <- "hgnc"
    familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData2, "_entrezID_family_TAD_DT.Rdata")))))
    familyDT$entrezID <- as.character(familyDT$entrezID)
    
    tadpos_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt")
    tadpos_dt <- read.delim(tadpos_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo","region", "start", "end"))
    stopifnot(is.numeric(tadpos_dt$start))
    stopifnot(is.numeric(tadpos_dt$end))
    tadpos_dt$region <- as.character(tadpos_dt$region)
    stopifnot(tadpos_dt$end > tadpos_dt$start)
    tadpos_dt <- tadpos_dt[order(tadpos_dt$chromo, tadpos_dt$start),]
    stopifnot(!duplicated(tadpos_dt$region))
    tadpos_dt$tad_rank <- 1:nrow(tadpos_dt)
    tad2rank <- setNames(tadpos_dt$tad_rank, tadpos_dt$region)      
                
    g2t_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt")
    g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    tad_g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    stopifnot(!duplicated(g2t_dt$entrezID))
    entrezIDchromo <- setNames(g2t_dt$chromo, g2t_dt$entrezID)
    
    # entrezIDstarts <- setNames(g2t_dt$start, g2t_dt$entrezID) # because of the strand -> done if gff_dt
    stopifnot(g2t_dt$entrezID %in% names(all_entrezIDstarts))
    entrezIDstarts <- all_entrezIDstarts[names(all_entrezIDstarts) %in% g2t_dt$entrezID]
    stopifnot(is.numeric(entrezIDstarts))
    entrezIDregion <- setNames(g2t_dt$region, g2t_dt$entrezID)
    stopifnot(entrezIDregion %in% names(tad2rank))
    
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
    
    # all_fams=all_fams[758]
    
    ##>> iterate here over families
    i_fam=1
    all_fam_results <- foreach(i_fam = 1:length(all_fams)) %dopar% {
      # all_fam_results <- foreach(i_fam = 1) %dopar% {
      
      fam <- all_fams[i_fam]
      
      cat(paste0("... ", hicds, " - start fam: ", i_fam, "/", length(all_fams), "\n"))
      
      fam_dt <- sameFamSameTAD_dt[sameFamSameTAD_dt$family == fam,]
      
      # sameTAD_fam_dt <- fam_dt[!is.na(fam_dt$region),]
      
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
      
      ## should be done chromo by chromo !!!!
      chromo_nbrBD_results <- foreach(chromo = names(true_countChromos)) %do% {
        
        chromo_tadpos_dt <- tadpos_dt[tadpos_dt$chromo == chromo,]
        chromo_size <- max(chromo_tadpos_dt$end)
        
        # v2: change here 19.10.2020 -> sample from all genes
        # all_chromo_genes <- names(family_entrezIDchromo)[family_entrezIDchromo == chromo]
        all_chromo_genes <- names(tad_entrezIDchromo)[tad_entrezIDchromo == chromo]
        true_chromo_genes <- true_genes[true_genes %in% all_chromo_genes]
        
        ### RETRIEVE THE CONTIGS AND THE GENES THEY SPAN
        true_ranks <- sort(entrezID_geneRanks[true_chromo_genes])
        stopifnot(length(true_ranks) == length(true_chromo_genes))
        cont_gene_ranks <- diff(true_ranks)
        stopifnot(cont_gene_ranks > 0)
        rle_cont_genes <- rle(cont_gene_ranks)
        contig_tosample <- rle_cont_genes$lengths[rle_cont_genes$values==1]+1  # this works even if empty : if c() + 1 = numeric(0)
        
        if(length(contig_tosample)== 0) return(NULL)
        
        contig_idx <- which(rle_cont_genes$values == 1)
        
        contig_ends <- cumsum(rle_cont_genes$lengths)[contig_idx] + 1
        # in fact it is (true-end - (lengths+1) + 1) => true-end - lengths
        contig_starts <- contig_ends -  rle_cont_genes$lengths[rle_cont_genes$values == 1]
        
        nContigs <- length(contig_ends)
        stopifnot(nContigs == length(contig_tosample))
        
        stopifnot(nContigs > 0)
        
        # for each contig:
        # -> retrieve the tad rank of the first gene
        # -> retrieve the tad rank of the last gene
        # count # of boundaries
        
        
        contig_results <- foreach(i = 1:nContigs) %do% {
          
          nGenes <- length(contig_starts[i]:contig_ends[i])
          stopifnot(nGenes == (rle_cont_genes$lengths[rle_cont_genes$values==1]+1)[i])
          
          stopifnot(diff(true_ranks[contig_starts[i]:contig_ends[i]]) == 1)
          gene_contig_first <-  names(true_ranks)[contig_starts[i]]
          gene_contig_last <-  names(true_ranks)[contig_ends[i]]
          
          stopifnot(gene_contig_first %in% names(entrezIDstarts))
          stopifnot(gene_contig_last %in% names(entrezIDstarts))
          stopifnot(gene_contig_first %in% names(entrezIDregion))
          stopifnot(gene_contig_last %in% names(entrezIDregion))
          
          
          # tad_firstGene <- entrezIDregion[gene_contig_first]
          # tad_lastGene <- entrezIDregion[gene_contig_last]
          # 
          # stopifnot(tad_firstGene %in% names(tad2rank))
          # stopifnot(tad_lastGene %in% names(tad2rank))
          
          # obs_nbrBD <- as.numeric(tad2rank[tad_lastGene]) -as.numeric(tad2rank[tad_firstGene])
          # stopifnot(is.numeric(obs_nbrBD))
          # stopifnot(obs_nbrBD >= 0)
          
          # need to do in this way to have the same way of assigning genes to tad
          
          startpos_firstGene <- as.numeric(entrezIDstarts[gene_contig_first])
          startpos_lastGene <- as.numeric(entrezIDstarts[gene_contig_last])
          contig_dist <- startpos_lastGene- startpos_firstGene + 1 # start=1, end=40, dist=40
          stopifnot(is.numeric(contig_dist))
          stopifnot(contig_dist > 0)
          
          obs_start <- startpos_firstGene
          obs_end <- startpos_lastGene
          
          obs_tad_firstGene <- chromo_tadpos_dt$region[chromo_tadpos_dt$start <= obs_start &
                                                            chromo_tadpos_dt$end >= obs_start]
          stopifnot(length(obs_tad_firstGene) == 1)
          
          obs_tad_lastGene <- chromo_tadpos_dt$region[chromo_tadpos_dt$start <= obs_end &
                                                           chromo_tadpos_dt$end >= obs_end]
          stopifnot(length(obs_tad_lastGene) == 1)
          
          obs_nbrBD <- as.numeric(tad2rank[obs_tad_lastGene]) - as.numeric(tad2rank[obs_tad_firstGene])
          stopifnot(is.numeric(obs_nbrBD))
          stopifnot(obs_nbrBD >= 0)
          
          
          
          
          
          contig_max_pos <- chromo_size - contig_dist+1 # if size=100, maxPos = 61
          
          random_nbrBD <- foreach(i = 1:nRandom, .combine='c') %do% {
            random_start <- sample(1:contig_max_pos, size=1)
            random_end <- random_start + contig_dist - 1 # if random start 61, random_end = 100
            stopifnot(random_end <=  chromo_size)
            # now retrieve the region of the random_start
            random_tad_firstGene <- chromo_tadpos_dt$region[chromo_tadpos_dt$start <= random_start &
                                                              chromo_tadpos_dt$end >= random_start]
            stopifnot(length(random_tad_firstGene) == 1)
            
            random_tad_lastGene <- chromo_tadpos_dt$region[chromo_tadpos_dt$start <= random_end &
                                                              chromo_tadpos_dt$end >= random_end]
            stopifnot(length(random_tad_lastGene) == 1)
            
            rd_nbrBD <- as.numeric(tad2rank[random_tad_lastGene]) - as.numeric(tad2rank[random_tad_firstGene])
            stopifnot(is.numeric(rd_nbrBD))
            stopifnot(rd_nbrBD >= 0)
            rd_nbrBD
          } #end-foreach iterating over number of permutations

          list(
            contig_dist=contig_dist,
            nGenes=nGenes,
            obs_nbrBD = obs_nbrBD,
            random_nbrBD = random_nbrBD
          )
          
        } # end-foreach iterating over the contigs
        names(contig_results) <- paste0("contig", 1:nContigs)
        contig_results
          
      } # end-foreach iterating over the chromos
      names(chromo_nbrBD_results) <- names(true_countChromos)
      chromo_nbrBD_results    
                    
    } # end-foreach iterating over the families
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
all_ds_results <- get(load("SAMEFAMCLUSTERNBRBD/all_ds_results.Rdata"))



# 1st level 
nGenes_dt <- lapply(all_ds_results, function(sub1) 
              lapply(sub1, function(sub2)
                lapply(sub2, function(sub3)
                  lapply(sub3, function(x)x[["nGenes"]]))))
nGenes_dt <- data.frame(unlist(nGenes_dt))
nGenes_dt$id <- rownames(nGenes_dt)
colnames(nGenes_dt) <- c("nGenes", "id")


obsBD_dt <- lapply(all_ds_results, function(sub1) 
  lapply(sub1, function(sub2)
    lapply(sub2, function(sub3)
      lapply(sub3, function(x)x[["obs_nbrBD"]]))))
obsBD_dt <- data.frame(unlist(obsBD_dt))
obsBD_dt$id <- rownames(obsBD_dt)
colnames(obsBD_dt) <- c("obs_nbrBD", "id")


meanrdBD_dt <- lapply(all_ds_results, function(sub1) 
  lapply(sub1, function(sub2)
    lapply(sub2, function(sub3)
      lapply(sub3, function(x)mean(x[["random_nbrBD"]])))))
meanrdBD_dt <- data.frame(unlist(meanrdBD_dt))
meanrdBD_dt$id <- rownames(meanrdBD_dt)
colnames(meanrdBD_dt) <- c("mean_random_nbrBD", "id")


merged_dt <- merge(meanrdBD_dt, merge(nGenes_dt, obsBD_dt, by="id", all=TRUE), by="id", all=TRUE)

plot(merged_dt$mean_random_nbrBD~merged_dt$obs_nbrBD)
curve(1*x, add=T)

count_obs_dt <- data.frame(
  type="observed",
  nbrBD= as.character(names(table(merged_dt$obs_nbrBD))),
  count= as.numeric(table(merged_dt$obs_nbrBD)),
  ratio= as.numeric(table(merged_dt$obs_nbrBD))/sum(table(merged_dt$obs_nbrBD)),
  stringsAsFactors = FALSE)

# round for the random !
rdd_rd_data <- round(merged_dt$mean_random_nbrBD)
count_rd_dt <- data.frame(
  type="mean_random",
  nbrBD= as.character(names(table(rdd_rd_data))),
  count= as.numeric(table(rdd_rd_data)),
  ratio= as.numeric(table(rdd_rd_data))/sum(table(rdd_rd_data)),
  stringsAsFactors = FALSE)

plot_dt <- rbind(count_obs_dt, count_rd_dt)

require(ggpubr)
require(ggsci)
plotTit <- ""
subTit <- ""

ggbarplot(plot_dt, x = "type", y="ratio", 
          fill = "nbrBD", 
          xlab="",
          ylab = "") + 
  scale_fill_nejm() + 
  labs(fill="") + 
  ggtitle(plotTit, subtitle=subTit)+
  theme(
    plot.title = element_text(size=16, face = "bold", hjust=0.5),
    plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
  )

outFile <- file.path(outFolder, paste0("nbrUniqueTADsByFam_log10_by_cat_barplot.", plotType))
ggsave(p1_nbr, filename=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile,  "\n"))





###################################################################### toy data checking contig retrieval
# #A)
# true_ranks <- c(
#   "g1" = 5,
#   "g2" = 6,
#   "g3" = 7,
#   "g4" = 10,
#   "g5" = 12,
#   "g6" = 13,
#   "g7" = 15
# )
# #B)
# true_ranks <- c(
#   "g0" = 1,
#   "g1" = 5,
#   "g2" = 6,
#   "g3" = 7,
#   "g4" = 10,
#   "g5" = 12,
#   "g6" = 13,
#   "g7" = 15
# )
# #C)
# true_ranks <- c(
#   "g0" = 1,
#   "g1" = 5,
#   "g2" = 6,
#   "g3" = 7,
#   "g4" = 10,
#   "g5" = 12,
#   "g6" = 13,
#   "g7" = 16,
#   "g8" = 17,
#   "g9" = 18
# )
# #D)
# true_ranks <- c(
#   "g0" = 100,
#   "g1" = 101,
#   "g2" = 504,
#   "g3" = 505,
#   "g4" = 506,
#   "g5" = 1218,
#   "g6" = 10000,
#   "g7" = 10001,
#   "g8" = 200000
# )
# 
# cont_gene_ranks <- diff(true_ranks)
# stopifnot(cont_gene_ranks > 0)
# rle_cont_genes <- rle(cont_gene_ranks)
# contig_tosample <- rle_cont_genes$lengths[rle_cont_genes$values==1]+1
# 
# contig_idx <- which(rle_cont_genes$values == 1)
# 
# contig_ends <- cumsum(rle_cont_genes$lengths)[contig_idx] + 1
# # in fact it is (true-end - (lengths+1) + 1) => true-end - lengths
# contig_starts <- contig_ends -  rle_cont_genes$lengths[rle_cont_genes$values == 1]
# 
# nContigs <- length(contig_ends)
# stopifnot(nContigs == length(contig_tosample))
# 
# all_contigs <- lapply(1:nContigs, function(i) {
#   stopifnot(diff(true_ranks[contig_starts[i]:contig_ends[i]]) == 1)
#   gene_contig_start <-  names(true_ranks)[contig_starts[i]]
#   gene_contig_end <-  names(true_ranks)[contig_ends[i]]
#   c(gene_contig_start, gene_contig_end)
# })
# all_contigs
# # A) [g1,g3], [g5,g6]
# # B) [g1,g3], [g5,g6]
# # C) [g1,g3], [g5, g6], [g7,g9]
# # D) [g0,g1], [g2, g4], [g6,g7]

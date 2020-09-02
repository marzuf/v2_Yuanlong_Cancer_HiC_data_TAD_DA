

# Rscript sameFamInSameTAD.R

# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity

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


plotType <- "svg"
myHeight <- 5
myWidth <- 7

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- file.path("SAMEFAMINSAMETAD")
dir.create(outFolder, recursive = TRUE)

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))


all_ds <- unlist(sapply(all_hicds, function(x) file.path(x, list.files(file.path(pipFolder, x)))), use.names = FALSE)

hicds = "Barutcu_MCF-10A_40kb"
exprds = "TCGAbrca_lum_bas"
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]

buildData <- TRUE

# all_ds=all_ds[1]

if(buildData) {
  
  
  all_ds_results <- foreach(ds = all_ds) %do%{
    
    hicds <- dirname(ds)
    exprds <- basename(ds)
    cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
    
    
    #       ### => CHANGED FOR THE TISSUE DATA TO USE TISSUE SPECIFIC FAMILY FILES !!!
    #       # inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
          inFoldFamily <- file.path("PREP_GENE_FAMILIES_TAD_DATA", hicds)
          familyData2 <- "hgnc"
          familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData2, "_entrezID_family_TAD_DT.Rdata")))))
          familyDT$entrezID <- as.character(familyDT$entrezID)
    
    g2t_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt")
    g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(!duplicated(g2t_dt$entrezID))
    entrezIDchromo <- setNames(g2t_dt$chromo, g2t_dt$entrezID)

    sameTADfile <- file.path("CREATE_SAME_TAD_SORTNODUP", hicds, "all_TAD_pairs.Rdata")
    stopifnot(file.exists(sameTADfile))
    sameTAD_dt <- get(load(sameTADfile))
    
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
    
    sameFamSameTAD_dt <- merge(sameFam_dt, sameTAD_dt, all.x=TRUE, all.y=FALSE, by=c("gene1", "gene2"))
    
    all_genes <- familyDT$entrezID
    stopifnot(all_genes %in% names(entrezIDchromo))
    family_entrezIDchromo <- entrezIDchromo[names(entrezIDchromo) %in% all_genes]
    
    
    all_fams <- unique(sameFamSameTAD_dt$family)
    
    ##>> iterate here over families
    
    all_fam_results <- foreach(i_fam = 1:length(all_fams)) %dopar% {
      
      fam <- all_fams[i_fam]
      
      cat(paste0("... ", hicds, " - ", exprds, " - start fam: ", i_fam, "/", length(all_fams), "\n"))
      
      fam_dt <- sameFamSameTAD_dt[sameFamSameTAD_dt$family == fam,]
      
      sameTAD_fam_dt <- fam_dt[!is.na(fam_dt$region),]
      
      nbrGenes <- length(unique(c(fam_dt$gene1, fam_dt$gene2)))
      
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
      
      random_results_dt <- foreach(i = 1:nRandom, .combine='cbind') %do% {
        
        sample_genes <- foreach(chromo = names(true_countChromos), .combine='c') %do% {
          sample(names(family_entrezIDchromo)[family_entrezIDchromo == chromo], size=true_countChromos[chromo], replace=FALSE)
        }
        stopifnot(!duplicated(sample_genes))
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
        
        random_fam_net <- graph_from_data_frame(d=random_sameTAD_fam_dt[,c("gene1", "gene2")], directed=F) 
        random_nbrEdges <- gsize(random_fam_net)
        stopifnot(random_nbrEdges == nrow(random_sameTAD_fam_dt))
        
        random_nbrComponents <- components(random_fam_net)$no
        stopifnot(random_nbrComponents == length(unique(random_sameTAD_fam_dt$region)))
        
        c(random_nbrEdges=random_nbrEdges, random_nbrComponents=random_nbrComponents)
      } # end permut
      
      list(
        nbrGenes=nbrGenes,
        nbrEdges=nbrEdges,
        nbrComponents=nbrComponents,
        random_results_dt=random_results_dt
      )
    } # end families
    names(all_fam_results) <- all_fams
    all_fam_results
  } # end datasets
  names(all_ds_results) <- all_ds
  
  outFile <- file.path(outFolder, "all_ds_results.Rdata")
  save(all_ds_results, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_ds_results.Rdata")
  all_ds_results <- get(load(outFile))
}
    
    
    
    
#     
#     
#     
#       
#       dataset_pipDir <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds) # used to retrieve gene list
#       stopifnot(dir.exists(dataset_pipDir))
#       
#       ### => CHANGED FOR THE TISSUE DATA TO USE TISSUE SPECIFIC FAMILY FILES !!!
#       # inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
#       inFoldFamily <- file.path("PREP_GENE_FAMILIES_TAD_DATA", hicds)
#       familyData2 <- "hgnc"
#       familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData2, "_entrezID_family_TAD_DT.Rdata")))))
#       familyDT$entrezID <- as.character(familyDT$entrezID)
#       
#       # entrezID chromo   start     end      region
#       # 49    10276  chr10 5454514 5501019 chr10_TAD17
#       # 50    51806  chr10 5540658 5541533 chr10_TAD18
#       # hgnc_family
#       # 49 Pleckstrin homology domain containing|Rho guanine nucleotide exchange factors
#       # 50                                                     EF-hand domain containing
#       # hgnc_family_short
#       
#       # -> for each family
#       # -> subset entrez
#       # -> retrieve the ones in dist <= meantadsize
#       # -> check: meanCorr == average of coexpr from coexprdt
#       
#       cat(paste0("... load DIST data\t", distFile, "\t", Sys.time(), "\t"))
#       load(distFile)
#       cat(paste0(Sys.time(), "\n"))
#       head(all_dist_pairs)
#       nrow(all_dist_pairs)
#       all_dist_pairs$gene1 <- as.character(all_dist_pairs$gene1)
#       all_dist_pairs$gene2 <- as.character(all_dist_pairs$gene2)
#       # UPDATE 30.06.2018
#       stopifnot(all_dist_pairs$gene1 < all_dist_pairs$gene2)
#       
#       # cat(paste0("... load TAD data\t", sameTADfile, "\t", Sys.time(), "\t"))
#       # ### =>>> CHANGED HERE FOR OTHER TAD FILE !!!
#       # load(sameTADfile)
#       # cat(paste0(Sys.time(), "\n"))
#       # head(all_TAD_pairs)
#       # nrow(all_TAD_pairs)
#       # all_TAD_pairs$gene1 <- as.character(all_TAD_pairs$gene1)
#       # all_TAD_pairs$gene2 <- as.character(all_TAD_pairs$gene2)
#       # # UPDATE 30.06.2018
#       # stopifnot(all_TAD_pairs$gene1 < all_TAD_pairs$gene2)
#       
#       # cat(paste0("... load COEXPR data\t",coexprFile, "\t", Sys.time(), "\t"))
#       # load(coexprFile)
#       # cat(paste0(Sys.time(), "\n"))
#       # head(coexprDT)
#       # nrow(coexprDT)
#       # coexprDT$gene1 <- as.character(coexprDT$gene1)
#       # coexprDT$gene2 <- as.character(coexprDT$gene2)
#       # all_TAD_pairs$gene2
#       # # UPDATE 30.06.2018
#       # stopifnot(coexprDT$gene1 < coexprDT$gene2)
#       all_fams <- unique(familyDT[,familyVar])
#       # myfam = all_fams[1]
#       all_fams_dt <-  foreach(myfam = all_fams) %dopar% {
#         
#         cat(paste0(hicds, " - " , exprds, " - " , myfam, "\n"))
#         
#         myfam_entrez <- familyDT$entrezID[familyDT[,familyVar] == myfam]
#         
#         stopifnot(myfam_entrez %in% all_dist_pairs$gene1 |  myfam_entrez %in% all_dist_pairs$gene2)
#         # stopifnot(myfam_entrez %in% all_TAD_pairs$gene1 |  myfam_entrez %in% all_TAD_pairs$gene2)
#         
#         fam_dist_pairs <- all_dist_pairs[all_dist_pairs$gene1 %in% myfam_entrez &
#                                            all_dist_pairs$gene2 %in% myfam_entrez, ]
#         
#         
#         fam_edge_table_tmp <- fam_dist_pairs[fam_dist_pairs$dist <= maxGeneDist,c("gene1","gene2", "dist")]
#         fam_edge_table <- fam_edge_table_tmp
#         fam_edge_table$dist <- NULL
#         stopifnot(fam_edge_table$gene1 < fam_edge_table$gene2)
#         
#         fam_net <- graph_from_data_frame(d=fam_edge_table, directed=F) 
#         # plot it
#         # plot(fam_net)
#         fam_maxCliques <-  max_cliques(fam_net)
#         cl_fam_to_keep <- which(lengths(fam_maxCliques) >= minCliqueSize)
#         if(length(cl_fam_to_keep) == 0) return(NULL)
#         fam_maxCliques <- fam_maxCliques[cl_fam_to_keep]
#         
#         fam_maxCliques_dt <- do.call(rbind, lapply(1:length(fam_maxCliques), function(i) {
#           # cl_genes <- names(fam_maxCliques[[i]])
#           # gene_cmbs_dt <- as.data.frame(t(combn(cl_genes, m=2)))
#           # gene_cmbs_dt$V1 <- as.character(gene_cmbs_dt$V1)
#           # gene_cmbs_dt$V2 <- as.character(gene_cmbs_dt$V2)
#           # colnames(gene_cmbs_dt) <- c("gene1", "gene2")
#           # gene_cmbs_dt$clique <- paste0(myfam, "_cl", i) 
#           # gene_cmbs_dt
#           data.frame(
#             entrezID = names(fam_maxCliques[[i]]),
#             clique = paste0(myfam, "_cl", i) ,
#             stringsAsFactors = FALSE
#           )
#         }))
#         stopifnot(!is.na(fam_maxCliques_dt))
#         stopifnot(table(fam_maxCliques_dt$clique) > 1)
#         
#         fam_maxCliques_dt$entrezID <- as.character(fam_maxCliques_dt$entrezID)
#         fam_maxCliques_dt$clique <- as.character(fam_maxCliques_dt$clique)
#         stopifnot(table(fam_maxCliques_dt$clique) >= minCliqueSize)
#         
#         # retrieve the maxdist
#         cl_tmp_dt <- data.frame(do.call(rbind, by(fam_maxCliques_dt, fam_maxCliques_dt$clique,function(x) t(combn(x$entrezID,2)))))
#         colnames(cl_tmp_dt) <- c("gene1", "gene2")
#         
#         stopifnot(cl_tmp_dt$gene1 %in% fam_edge_table_tmp$gene1 | cl_tmp_dt$gene1 %in% fam_edge_table_tmp$gene2)
#         stopifnot(cl_tmp_dt$gene2 %in% fam_edge_table_tmp$gene1 | cl_tmp_dt$gene2 %in% fam_edge_table_tmp$gene2)
#         
#         tmp_dt <- fam_edge_table_tmp[
#           (fam_edge_table_tmp$gene1 %in% cl_tmp_dt$gene1 |  fam_edge_table_tmp$gene1 %in% cl_tmp_dt$gene2) &
#             (fam_edge_table_tmp$gene2 %in% cl_tmp_dt$gene1 |  fam_edge_table_tmp$gene2 %in% cl_tmp_dt$gene2),
#         ]
#         maxCliqueDist <- as.numeric(max(tmp_dt$dist))
#         stopifnot(!is.na(maxCliqueDist))
#         
#         list(
#           maxPairDist = maxGeneDist,
#           maxCliqueDist=maxCliqueDist,
#           fam_cl_dt = fam_maxCliques_dt
#         )
#         
#       }
#       
#       outFile <- file.path(outFolder, hicds, exprds, "all_fams_dt.Rdata")
#       dir.create(dirname(outFile), recursive = TRUE)
#       save(all_fams_dt, file = outFile,version=2 )
#       cat(paste0("... written: ", outFile, "\n"))
#       
#     }
#   }
#   
#   
# } 
# 
# 
# maxCliqueDist_dt <- foreach(hicds = all_hicds, .combine='rbind') %do%{
#   cat(paste0("... start: ", hicds, "\n"))
#   exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
#     cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
#     
#     # retrieve file 
#     famMod_file <- file.path(outFolder, hicds, exprds, "all_fams_dt.Rdata")
#     stopifnot(file.exists(famMod_file))
#     fam_data <- get(load(famMod_file))
#     
#     all_maxCliqueDist <- as.numeric(unlist(lapply(fam_data, function(x) x[["maxCliqueDist"]])))
#     
#     data.frame(
#       hicds=hicds,
#       exprds=exprds,
#       maxCliqueDist=all_maxCliqueDist,
#       stringsAsFactors = FALSE
#     )
#   }
#   exprds_dt
# }
# 
# maxPairDist_dt <- foreach(hicds = all_hicds, .combine='rbind') %do%{
#   cat(paste0("... start: ", hicds, "\n"))
#   exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
#     cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
#     
#     # retrieve file 
#     famMod_file <- file.path(outFolder, hicds, exprds, "all_fams_dt.Rdata")
#     stopifnot(file.exists(famMod_file))
#     fam_data <- get(load(famMod_file))
#     
#     all_maxPairDist <- unique(as.numeric(unlist(lapply(fam_data, function(x) x[["maxPairDist"]]))))
#     
#     data.frame(
#       hicds=hicds,
#       exprds=exprds,
#       maxPairDist=all_maxPairDist,
#       stringsAsFactors = FALSE
#     )
#   }
# }
# 
# 
# nDS <- length(unique(file.path(maxPairDist_dt$hicds, maxPairDist_dt$exprds)))
# stopifnot(nDS == length(unique(file.path(maxCliqueDist_dt$hicds, maxCliqueDist_dt$exprds))))
# 
# 
# outFile <- file.path(outFolder,  paste0("allDS_distCliquesGenes_density.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# plot_multiDens(
#   list(maxPairDist = (maxPairDist_dt$maxPairDist),
#        maxCliqueDist = (maxCliqueDist_dt$maxCliqueDist)),
#   plotTit = paste0("all datasets -  n=", nDS)
# )
# mtext(side=3, text = paste0("minCliqueSize=", minCliqueSize), font=3)
# foo <- dev.off()
# cat(paste0("... written: ", outFile,  "\n"))
# 
# 
# outFile <- file.path(outFolder,  paste0("allDS_distCliquesGenes_density_log10.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# plot_multiDens(
#   list(maxPairDist_log10 = log10(maxPairDist_dt$maxPairDist),
#        maxCliqueDist_log10 = log10(maxCliqueDist_dt$maxCliqueDist)),
#   plotTit = paste0("all datasets -  n=", nDS)
# )
# mtext(side=3, text = paste0("minCliqueSize=", minCliqueSize), font=3)
# foo <- dev.off()
# cat(paste0("... written: ", outFile,  "\n"))

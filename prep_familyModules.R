

# Rscript prep_familyModules.R 

# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity

require(doMC)
require(foreach)
registerDoMC(40)
require(reshape2)
require(igraph)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

runFolder <- "."
pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")
familyVar <- "hgnc_family_short"
minCmpntSize <- 3

nMaxSize <- 2

outFolder <- file.path("PREP_FAMILYMODULES", nMaxSize)

plotType <- "svg"
myHeight <- 5
myWidth <- 7

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

hicds = "Barutcu_MCF-10A_40kb"
exprds="TCGAbrca_lum_bas"

# all_hicds=all_hicds[2:length(all_hicds)]
# all_hicds=all_hicds[1]

buildData <- TRUE

if(buildData){
  
  FOO1 <- foreach(hicds = all_hicds) %do%{
    cat(paste0("... start: ", hicds, "\n"))
    FOO2 <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
      
      # fcc_file <- file.path(pipFolder, hicds, exprds, step8fcc_folder, "all_obs_prodSignedRatio.Rdata")
      # stopifnot(file.exists(fcc_file))  
      
      tad_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt")
      tad_dt <- read.delim(tad_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo","region", "start", "end"))
      
      onlyTAD_dt <- tad_dt[grep("_TAD", tad_dt$region),]
      
      onlyTAD_dt$midPos <- 0.5*(onlyTAD_dt$end + onlyTAD_dt$start)
      
      # retrieve the average TAD size
#      meanTADsize <- mean(onlyTAD_dt$midPos) # 14.05.2020 CORRECTED
      meanTADsize <- mean(onlyTAD_dt$end-onlyTAD_dt$start+1) # 14.05.2020 CORRECTED
      stopifnot(!is.na(meanTADsize))
      
      maxGeneDist <- nMaxSize*meanTADsize
      
      corMethod <- "pearson"
      familyData <- "hgnc_family_short"
      
      ### UPDATE 08.01.19 => USE TAD LIST SPECIFIC FILES = TISSUE SPECIFIC FILES
      distFile <- file.path("CREATE_DIST_SORTNODUP", hicds, "all_dist_pairs.Rdata")
      stopifnot(file.exists(distFile))
      
      # CHANGED 08.01.19 !!!
      # coexprFile <- file.path("CREATE_COEXPR_SORTNODUP", hicds, exprds,  corMethod, "coexprDT.Rdata")
      # stopifnot(file.exists(coexprFile))
      
      # sameTADfile <- file.path("CREATE_SAME_TAD_SORTNODUP", hicds, "all_TAD_pairs.Rdata")
      # stopifnot(file.exists(sameTADfile))
      
      # ADDED 08.01.19 to accommodate updated family file
      sameFamFolder <- file.path("CREATE_SAME_FAMILY_SORTNODUP", hicds)
      # checking the file comes after (iterating over family and family_short)
      stopifnot(dir.exists(sameFamFolder))
      sameFamFile <- file.path(sameFamFolder, paste0(familyData, "_all_family_pairs.Rdata")) # at least this one should exist !
      stopifnot(file.exists(sameFamFile))
      
      dataset_pipDir <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds) # used to retrieve gene list
      stopifnot(dir.exists(dataset_pipDir))
      
      ### => CHANGED FOR THE TISSUE DATA TO USE TISSUE SPECIFIC FAMILY FILES !!!
      # inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
      inFoldFamily <- file.path("PREP_GENE_FAMILIES_TAD_DATA", hicds)
      familyData2 <- "hgnc"
      familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData2, "_entrezID_family_TAD_DT.Rdata")))))
      familyDT$entrezID <- as.character(familyDT$entrezID)
      
      # entrezID chromo   start     end      region
      # 49    10276  chr10 5454514 5501019 chr10_TAD17
      # 50    51806  chr10 5540658 5541533 chr10_TAD18
      # hgnc_family
      # 49 Pleckstrin homology domain containing|Rho guanine nucleotide exchange factors
      # 50                                                     EF-hand domain containing
      # hgnc_family_short
      
      # -> for each family
      # -> subset entrez
      # -> retrieve the ones in dist <= meantadsize
      # -> check: meanCorr == average of coexpr from coexprdt
      
      cat(paste0("... load DIST data\t", distFile, "\t", Sys.time(), "\t"))
      load(distFile)
      cat(paste0(Sys.time(), "\n"))
      head(all_dist_pairs)
      nrow(all_dist_pairs)
      all_dist_pairs$gene1 <- as.character(all_dist_pairs$gene1)
      all_dist_pairs$gene2 <- as.character(all_dist_pairs$gene2)
      # UPDATE 30.06.2018
      stopifnot(all_dist_pairs$gene1 < all_dist_pairs$gene2)
      
      # cat(paste0("... load TAD data\t", sameTADfile, "\t", Sys.time(), "\t"))
      # ### =>>> CHANGED HERE FOR OTHER TAD FILE !!!
      # load(sameTADfile)
      # cat(paste0(Sys.time(), "\n"))
      # head(all_TAD_pairs)
      # nrow(all_TAD_pairs)
      # all_TAD_pairs$gene1 <- as.character(all_TAD_pairs$gene1)
      # all_TAD_pairs$gene2 <- as.character(all_TAD_pairs$gene2)
      # # UPDATE 30.06.2018
      # stopifnot(all_TAD_pairs$gene1 < all_TAD_pairs$gene2)
      
      # cat(paste0("... load COEXPR data\t",coexprFile, "\t", Sys.time(), "\t"))
      # load(coexprFile)
      # cat(paste0(Sys.time(), "\n"))
      # head(coexprDT)
      # nrow(coexprDT)
      # coexprDT$gene1 <- as.character(coexprDT$gene1)
      # coexprDT$gene2 <- as.character(coexprDT$gene2)
      # all_TAD_pairs$gene2
      # # UPDATE 30.06.2018
      # stopifnot(coexprDT$gene1 < coexprDT$gene2)
      
      
      all_fams <- unique(familyDT[,familyVar])
      myfam = all_fams[1]
      # all_fams_dt <-  foreach(myfam = all_fams, .combine='rbind') %dopar% {
      all_fams_dt <-  foreach(myfam = all_fams) %dopar% {
        
        cat(paste0(hicds, " - " , exprds, " - " , myfam, "\n"))
        
        myfam_entrez <- familyDT$entrezID[familyDT[,familyVar] == myfam]
        
        stopifnot(myfam_entrez %in% all_dist_pairs$gene1 |  myfam_entrez %in% all_dist_pairs$gene2)
        # stopifnot(myfam_entrez %in% all_TAD_pairs$gene1 |  myfam_entrez %in% all_TAD_pairs$gene2)
        
        fam_dist_pairs <- all_dist_pairs[all_dist_pairs$gene1 %in% myfam_entrez &
                                           all_dist_pairs$gene2 %in% myfam_entrez, ]
        
        
        
        
        fam_edge_table_tmp <- fam_dist_pairs[fam_dist_pairs$dist <= maxGeneDist,c("gene1","gene2", "dist")]
        fam_edge_table <- fam_edge_table_tmp
        fam_edge_table$dist <- NULL
        stopifnot(fam_edge_table$gene1 < fam_edge_table$gene2)
        
        fam_net <- graph_from_data_frame(d=fam_edge_table, directed=F) 
        # plot it
        # plot(fam_net)
        fam_components <-  components(fam_net)$membership
        fam_size <-  components(fam_net)$csize
        fam_to_keep <- which(fam_size >= minCmpntSize)
        
        if(length(fam_to_keep) == 0) return(NULL)
        
        # fam_components_dt <- data.frame(
        #   entrezID=as.character(names(fam_components)),
        #   famCpmnt = as.numeric(fam_components),
        #   stringsAsFactors = FALSE
        # )
        # 347688  10376   7846  27175  51174 113691  51807 112714   7277 203068 347733 
        # 1      2      2      3      3      4      4      5      6      7      7 
        # 79861  84790   7283 6453 48 113457  80086  51175   7280 
        # 1      2      3      4      5      6      7      7 
        # fam_edge_coexpr_dt <- merge(fam_edge_table, coexprDT, by=c("gene1", "gene2"), all.x=TRUE, all.y=FALSE)
        # stopifnot(!is.na(fam_edge_coexpr_dt))
        fam_edge_cmpnt_dt <- fam_edge_table
        fam_edge_cmpnt_dt$cpt1 <- fam_components[as.character(fam_edge_cmpnt_dt$gene1)]
        fam_edge_cmpnt_dt$cpt2 <- fam_components[as.character(fam_edge_cmpnt_dt$gene2)]
        stopifnot(fam_edge_cmpnt_dt$cpt1 == fam_edge_cmpnt_dt$cpt2)
        
        keep_fam_edge_cmpnt_dt <- fam_edge_cmpnt_dt[fam_edge_cmpnt_dt$cpt1 %in% fam_to_keep,]
        disc_fam_edge_cmpnt_dt <- fam_edge_cmpnt_dt[! fam_edge_cmpnt_dt$cpt1 %in% fam_to_keep,]
        
        stopifnot(table(disc_fam_edge_cmpnt_dt$cpt1) <= 1)
        stopifnot(table(keep_fam_edge_cmpnt_dt$cpt1) > 1)
        
        out_dt <- keep_fam_edge_cmpnt_dt  
        out_dt$cpt <- paste0(myfam, "_cpt", out_dt$cpt2)
        out_dt$cpt2 <- NULL
        out_dt$cpt1 <- NULL
        
        stopifnot(out_dt$gene1 %in% fam_edge_table_tmp$gene1 | out_dt$gene1 %in% fam_edge_table_tmp$gene2)
        stopifnot(out_dt$gene2 %in% fam_edge_table_tmp$gene2 | out_dt$gene2 %in% fam_edge_table_tmp$gene2)
        
        # retrieve the actual max dist inside the component
        tmp_dt <- fam_edge_table_tmp[
          (fam_edge_table_tmp$gene1 %in% out_dt$gene1 |  fam_edge_table_tmp$gene1 %in% out_dt$gene2) &
            (fam_edge_table_tmp$gene2 %in% out_dt$gene1 |  fam_edge_table_tmp$gene2 %in% out_dt$gene2),
          ]
        maxCptDist <- as.numeric(max(tmp_dt$dist))
        stopifnot(!is.na(maxCptDist))
        
        out_dt_m <- melt(out_dt, id="cpt")
        out_dt_m$variable <- NULL
        out_dt_m <- unique(out_dt_m)
        stopifnot(table(out_dt_m$cpt) >= minCmpntSize)
        rownames(out_dt_m) <- NULL
        out_dt_m$entrezID <- as.character(out_dt_m$value)
        out_dt_m$value <- NULL
        stopifnot(!duplicated(out_dt_m$entrezID))
        list(
          fam_cpt_dt = out_dt_m,
          maxPairDist = maxGeneDist,
          maxCptDist = maxCptDist
        )
        
      }
      
      outFile <- file.path(outFolder, hicds, exprds, "all_fams_dt.Rdata")
      dir.create(dirname(outFile), recursive = TRUE)
      save(all_fams_dt, file = outFile,version=2 )
      cat(paste0("... written: ", outFile, "\n"))
      
    }
  }
  
  
} 


maxCptDist_dt <- foreach(hicds = all_hicds, .combine='rbind') %do%{
  cat(paste0("... start: ", hicds, "\n"))
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
    
    # retrieve file 
    famMod_file <- file.path(outFolder, hicds, exprds, "all_fams_dt.Rdata")
    stopifnot(file.exists(famMod_file))
    fam_data <- get(load(famMod_file))
    
    all_maxCptDist <- as.numeric(unlist(lapply(fam_data, function(x) x[["maxCptDist"]])))
    
    data.frame(
      hicds=hicds,
      exprds=exprds,
      maxCptDist=all_maxCptDist,
      stringsAsFactors = FALSE
    )
  }
  exprds_dt
}

maxPairDist_dt <- foreach(hicds = all_hicds, .combine='rbind') %do%{
  cat(paste0("... start: ", hicds, "\n"))
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
    
    # retrieve file 
    famMod_file <- file.path(outFolder, hicds, exprds, "all_fams_dt.Rdata")
    stopifnot(file.exists(famMod_file))
    fam_data <- get(load(famMod_file))
    
    all_maxPairDist <- unique(as.numeric(unlist(lapply(fam_data, function(x) x[["maxPairDist"]]))))
    
    data.frame(
      hicds=hicds,
      exprds=exprds,
      maxPairDist=all_maxPairDist,
      stringsAsFactors = FALSE
    )
  }
}


nDS <- length(unique(file.path(maxPairDist_dt$hicds, maxPairDist_dt$exprds)))
stopifnot(nDS == length(unique(file.path(maxCptDist_dt$hicds, maxCptDist_dt$exprds))))


outFile <- file.path(outFolder,  paste0("allDS_distCliquesGenes_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(maxPairDist = (maxPairDist_dt$maxPairDist),
       maxCptDist = (maxCptDist_dt$maxCptDist)),
  plotTit = paste0("all datasets -  n=", nDS)
)
mtext(side=3, text = paste0("minCmpntSize=", minCmpntSize), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))


outFile <- file.path(outFolder,  paste0("allDS_distCliquesGenes_density_log10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(maxPairDist_log10 = log10(maxPairDist_dt$maxPairDist),
       maxCptDist_log10 = log10(maxCptDist_dt$maxCptDist)),
  plotTit = paste0("all datasets =", nDS)
)
mtext(side=3, text = paste0("minCmpntSize=", minCmpntSize), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))



# ALTERNATIVE - but then tricky to recover actual maxdist:
# fam_components_dt <- data.frame(
#   entrezID=as.character(names(fam_components)),
#   famCpmnt = as.numeric(fam_components),
#   # famCpmnt = paste0(myfam, "_cpt", as.numeric(fam_components)),
#   stringsAsFactors = FALSE
# )
# # 347688  10376   7846  27175  51174 113691  51807 112714   7277 203068 347733 
# # 1      2      2      3      3      4      4      5      6      7      7 
# # 79861  84790   7283 6453 48 113457  80086  51175   7280 
# # 1      2      3      4      5      6      7      7 
# 
# keep_fam_cpt_dt <- fam_components_dt[ as.character(fam_components_dt$famCpmnt) %in% as.character(fam_to_keep),]
# stopifnot(nrow(keep_fam_cpt_dt) > 0)
# stopifnot(table(keep_fam_cpt_dt$famCpmnt) > 1)
# 
# keep_fam_cpt_dt$cpt <- paste0(myfam, "_cpt", as.character(keep_fam_cpt_dt$famCpmnt))
# 
# out_dt <- keep_fam_edge_cmpnt_dt  
# stopifnot(table(out_dt$cpt) >= minCmpntSize)
# out_dt$entrezID <- as.character(out_dt$entrezID)
# stopifnot(!duplicated(out_dt$entrezID))
# 
# # retrieve maxGeneDist
# stopifnot(out_dt$en)
# fam_edge_table_tmp

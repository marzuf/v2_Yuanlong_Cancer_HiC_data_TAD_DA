

# Rscript prep_familyModules.R

# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity

require(doMC)
require(foreach)
registerDoMC(40)

runFolder <- "."
binSize <- 40*10^3


all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]

hicds = "Barutcu_MCF-10A_40kb"


all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
  cat(paste0("... start: ", hicds, "\n"))
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
    
    
    
    # fcc_file <- file.path(pipFolder, hicds, exprds, step8fcc_folder, "all_obs_prodSignedRatio.Rdata")
    # stopifnot(file.exists(fcc_file))  

  
  
    tad_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt")
    tad_dt <- read.delim(tad_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo","region", "start", "end"))
    
    onlyTAD_dt <- tad_dt[grep("_TAD", tad_dt$region),]
    
    
    onlyTAD_dt$midPos <- 0.5*(onlyTAD_dt$end + onlyTAD_dt$start)
    
    # retrieve the average TAD size
    meanTADsize <- mean(onlyTAD_dt$midPos)
    stopifnot(!is.na(meanTADsize))
    
    corMethod <- "pearson"
    familyData <- "hgnc_family_short"
    
    ### UPDATE 08.01.19 => USE TAD LIST SPECIFIC FILES = TISSUE SPECIFIC FILES
    distFile <- file.path("CREATE_DIST_SORTNODUP", hicds, "all_dist_pairs.Rdata")
    stopifnot(file.exists(distFile))
    
    # CHANGED 08.01.19 !!!
    coexprFile <- file.path("CREATE_COEXPR_SORTNODUP", hicds, exprds,  corMethod, "coexprDT.Rdata")
    stopifnot(file.exists(coexprFile))
    
    sameTADfile <- file.path("CREATE_SAME_TAD_SORTNODUP", hicds, "all_TAD_pairs.Rdata")
    stopifnot(file.exists(sameTADfile))
    
    # ADDED 08.01.19 to accommodate updated family file
    sameFamFolder <- file.path("CREATE_SAME_FAMILY_SORTNODUP", hicds)
    # checking the file comes after (iterating over family and family_short)
    stopifnot(dir.exists(sameFamFolder))
    sameFamFile <- file.path(sameFamFolder, paste0(familyData, "_family_all_family_pairs.Rdata")) # at least this one should exist !
    stopifnot(file.exists(sameFamFile))
    
    dataset_pipDir <- file.path("PIPELINE", "OUTPUT_FOLDER", curr_TADlist, curr_dataset) # used to retrieve gene list
    stopifnot(dir.exists(dataset_pipDir))
    
    ### => CHANGED FOR THE TISSUE DATA TO USE TISSUE SPECIFIC FAMILY FILES !!!
    # inFoldFamily <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller, "/", "PREP_GENE_FAMILIES_TAD_DATA"))
    inFoldFamily <- file.path("PREP_GENE_FAMILIES_TAD_DATA", hicds)
    familyDT <- eval(parse(text = load(file.path(inFoldFamily, paste0(familyData, "_entrezID_family_TAD_DT.Rdata")))))
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
    
    txt <- paste0("... distFile = ",  distFile, "\n")
    
    
    
    
    
    cat(paste0("... load DIST data\t", distFile, "\t", Sys.time(), "\t"))
    load(distFile)
    cat(paste0(Sys.time(), "\n"))
    head(all_dist_pairs)
    nrow(all_dist_pairs)
    all_dist_pairs$gene1 <- as.character(all_dist_pairs$gene1)
    all_dist_pairs$gene2 <- as.character(all_dist_pairs$gene2)
    # UPDATE 30.06.2018
    stopifnot(all_dist_pairs$gene1 < all_dist_pairs$gene2)
    
    cat(paste0("... load TAD data\t", sameTADfile, "\t", Sys.time(), "\t"))
    ### =>>> CHANGED HERE FOR OTHER TAD FILE !!!
    load(sameTADfile)
    cat(paste0(Sys.time(), "\n"))
    head(all_TAD_pairs)
    nrow(all_TAD_pairs)
    all_TAD_pairs$gene1 <- as.character(all_TAD_pairs$gene1)
    all_TAD_pairs$gene2 <- as.character(all_TAD_pairs$gene2)
    # UPDATE 30.06.2018
    stopifnot(all_TAD_pairs$gene1 < all_TAD_pairs$gene2)
    
    cat(paste0("... load COEXPR data\t",coexprFile, "\t", Sys.time(), "\t"))
    load(coexprFile)
    cat(paste0(Sys.time(), "\n"))
    head(coexprDT)
    nrow(coexprDT)
    coexprDT$gene1 <- as.character(coexprDT$gene1)
    coexprDT$gene2 <- as.character(coexprDT$gene2)
    all_TAD_pairs$gene2
    # UPDATE 30.06.2018
    stopifnot(coexprDT$gene1 < coexprDT$gene2)
    
    
    
    
    
    
    foo <- foreach(chr = unique(tad_dt$chromo)) %dopar% {
      
      sub_dt <- onlyTAD_dt[onlyTAD_dt$chromo == chr,c("chromo", "start", "end", "midPos_bin")]
      
      stopifnot(sub_dt$end > sub_dt$start)
      stopifnot(sub_dt$end %% binSize == 0)
      stopifnot(sub_dt$start %% binSize == 1)
      
      new_tad_dt <- data.frame(
        chromo = chr,
        start = sub_dt$midPos_bin[1:(nrow(sub_dt)-1)]+1,
        end = sub_dt$midPos_bin[2:nrow(sub_dt)],
        stringsAsFactors = FALSE
      )
      
      stopifnot(new_tad_dt$end > new_tad_dt$start)
      stopifnot(new_tad_dt$end %% binSize == 0)
      stopifnot(new_tad_dt$start %% binSize == 1)
      
      outFile <- file.path(outFolder, paste0(rd_hicds, "_", chr, "_YL_", binSize/1000, "kb_final_domains.txt"))
      write.table(new_tad_dt, file = outFile, sep="\t", col.names=F, row.names=F, quote=F, append=F )
      cat(paste0("... written: ", outFile, "\n"))
      
      
    }
    
  }
  
  # call assign_genes
  mycmd <- paste0("./3_assign_genes.sh ", rd_hicds)
  cat(paste0("> ", mycmd, "\n"))
  system(mycmd)
}






# Rscript sameFamInSameTAD_regionOcc.R

# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity

set.seed(20200903) # NB forgot the first time I launched it

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

outFolder <- file.path("SAMEFAMINSAMETAD_REGIONOCC")
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

ds=all_ds[1]

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
    
    tad_g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    stopifnot(familyDT$entrezID %in% tad_g2t_dt$entrezID)
    entrezIDchromo <- setNames(tad_g2t_dt$chromo, tad_g2t_dt$entrezID)
    entrezIDregion <- setNames(tad_g2t_dt$region, tad_g2t_dt$entrezID)
    stopifnot(!duplicated(familyDT$entrezID))
    
    familyDT$region <- entrezIDregion[familyDT$entrezID]
    stopifnot(!is.na(familyDT$region))
    
    all_genes <- familyDT$entrezID
    stopifnot(all_genes %in% names(entrezIDchromo))
    family_entrezIDchromo <- entrezIDchromo[names(entrezIDchromo) %in% all_genes]
    
    
    all_fams <- unique(familyDT[,familyData])
    
    i_fam=1
    all_fam_results <- foreach(i_fam = 1:length(all_fams)) %dopar% {
      
      
      fam <- all_fams[i_fam]
      
      cat(paste0("... ", hicds, " - ", exprds, " - start fam: ", i_fam, "/", length(all_fams), "\n"))
      
      sub_dt <- familyDT[familyDT[,familyData] == fam,]
      
      nbrUniqueRegions <- length(unique(sub_dt$region))
      nbrGenes <- length(unique(sub_dt$entrezID))
      stopifnot(nrow(sub_dt) == nbrGenes)
      regionOcc <- c(table(sub_dt$region))
      
      
      
      true_genes <- sub_dt$entrezID
      stopifnot(true_genes %in% names(entrezIDchromo))
      true_chromos <- entrezIDchromo[true_genes]
      true_countChromos <- table(true_chromos)
      
      random_results_data <- foreach(i = 1:nRandom) %do% {
        
        sample_genes <- foreach(chromo = names(true_countChromos), .combine='c') %do% {
          sample(names(family_entrezIDchromo)[family_entrezIDchromo == chromo], size=true_countChromos[chromo], replace=FALSE)
        }
        stopifnot(!duplicated(sample_genes))
        stopifnot(length(sample_genes) == length(true_genes))
        
        stopifnot(sample_genes %in% familyDT$entrezID)
        random_sub_dt <- familyDT[familyDT$entrezID %in% sample_genes,]
        
        random_nbrUniqueRegions <- length(unique(random_sub_dt$region))
        random_nbrGenes <- length(unique(random_sub_dt$entrezID))
        stopifnot(nrow(random_sub_dt) == random_nbrGenes)
        random_regionOcc <- c(table(random_sub_dt$region))
        

        list(
          random_nbrUniqueRegions=random_nbrUniqueRegions,
          random_nbrGenes=random_nbrGenes,
          random_ratioUniqueRegions=random_nbrUniqueRegions/random_nbrGenes,
          random_regionOcc=random_regionOcc
        )
      } # end permut
      
      list(
        nbrUniqueRegions=nbrUniqueRegions,
        nbrGenes=nbrGenes,
        ratioUniqueRegions=nbrUniqueRegions/nbrGenes,
        regionOcc=regionOcc,
        random_results_data=random_results_data
      )
      
    } # end over fams
    
    names(all_fam_results) <- all_fams
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", "all_fam_results.Rdata"))
    save(all_fam_results, file=outFile, version=2)
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
    
all_ds_results = get(load("SAMEFAMINSAMETAD_REGIONOCCS/all_ds_results.Rdata"))

all_obs_ratioUniqueRegions <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["ratioUniqueRegions"]])))
stopifnot(all_obs_ratioUniqueRegions >= 0 & all_obs_ratioUniqueRegions <= 1)
all_random_ratioUniqueRegions <- unlist(lapply(all_ds_results, function(dslist) lapply(dslist, function(famlist) lapply(famlist[["random_results_data"]], function(x) x[["random_ratioUniqueRegions"]]))))
stopifnot(all_random_ratioUniqueRegions >= 0 & all_random_ratioUniqueRegions <= 1)

plot_multiDens(
  list(
    all_obs_ratioUniqueRegions=all_obs_ratioUniqueRegions,
    all_random_ratioUniqueRegions=all_random_ratioUniqueRegions
  )
)

# ratio of TADs with more than one gene
all_obs_ratioRegionsMultipleGenes <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) {
  mean(x[["regionOcc"]] > 1)
}))) 

all_random_ratioRegionsMultipleGenes <- unlist(lapply(all_ds_results, function(dslist) lapply(dslist, function(famlist) lapply(famlist[["random_results_data"]], function(x) {
  mean(x[["random_regionOcc"]] > 1)
}))))

plot_multiDens(
  list(
    all_obs_ratioRegionsMultipleGenes=all_obs_ratioRegionsMultipleGenes,
    all_random_ratioRegionsMultipleGenes=all_random_ratioRegionsMultipleGenes
  )
)

plot_multiDens(
  list(
    all_obs_ratioRegionsMultipleGenes=all_obs_ratioRegionsMultipleGenes[all_obs_ratioRegionsMultipleGenes>0],
    all_random_ratioRegionsMultipleGenes=all_random_ratioRegionsMultipleGenes[all_random_ratioRegionsMultipleGenes>0]
  )
)


#  barplot # with 1 unique Reg, 2-3 reg, 4-5, 5
all_obs_ratioUniqueRegions <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["ratioUniqueRegions"]])))

dslist=all_ds_results[[1]]
famlist=dslist[[1]]

all_cat_nbrUniqueTADs <- do.call(rbind, lapply(all_ds_results, function(dslist) do.call(rbind, lapply(dslist, function(famlist) {
  obs_uniqueTADs <- famlist[["nbrUniqueRegions"]]
  rd_uniqueTADs <- unlist(lapply(famlist[["random_results_data"]], function(x) x[["random_nbrUniqueRegions"]]))
  
  obs_cat1 <- sum(obs_uniqueTADs == 1)
  obs_cat2 <- sum(obs_uniqueTADs == 2 | obs_uniqueTADs == 3 )
  obs_cat3 <- sum(obs_uniqueTADs == 4 | obs_uniqueTADs == 5 )
  obs_cat4 <- sum(obs_uniqueTADs > 5 )
  
  rd_cat1 <- sum(rd_uniqueTADs == 1)
  rd_cat2 <- sum(rd_uniqueTADs == 2 | rd_uniqueTADs == 3 )
  rd_cat3 <- sum(rd_uniqueTADs == 4 | rd_uniqueTADs == 5 )
  rd_cat4 <- sum(rd_uniqueTADs > 5 )
  
  data.frame(
    obs_cat1=obs_cat1,
    obs_cat2=obs_cat2,
    obs_cat3=obs_cat3,
    obs_cat4=obs_cat4,
    rd_cat1=rd_cat1,
    rd_cat2=rd_cat2,
    rd_cat3=rd_cat3,
    rd_cat4=rd_cat4,
    stringsAsFactors = FALSE
  )
}))) )
rownames(all_cat_nbrUniqueTADs) <- NULL



all_cat_nbrUniqueTADs2 <- do.call(rbind, lapply(all_ds_results, function(dslist) {
  
  tmp_res <- lapply(dslist, function(famlist) {
  obs_uniqueTADs <- famlist[["nbrUniqueRegions"]]
  rd_uniqueTADs <- unlist(lapply(famlist[["random_results_data"]], function(x) x[["random_nbrUniqueRegions"]]))
    list(obs_uniqueTADs=obs_uniqueTADs, rd_uniqueTADs=rd_uniqueTADs)
  })
  obs_uniqueTADs <- unlist(lapply(tmp_res, function(x) x[["obs_uniqueTADs"]]))
  rd_uniqueTADs <- unlist(lapply(tmp_res, function(x) x [["rd_uniqueTADs"]]))
  
  obs_cat1 <- sum(obs_uniqueTADs == 1)
  obs_cat2 <- sum(obs_uniqueTADs == 2 | obs_uniqueTADs == 3 )
  obs_cat3 <- sum(obs_uniqueTADs == 4 | obs_uniqueTADs == 5 )
  obs_cat4 <- sum(obs_uniqueTADs > 5 )
  
  rd_cat1 <- sum(rd_uniqueTADs == 1)
  rd_cat2 <- sum(rd_uniqueTADs == 2 | rd_uniqueTADs == 3 )
  rd_cat3 <- sum(rd_uniqueTADs == 4 | rd_uniqueTADs == 5 )
  rd_cat4 <- sum(rd_uniqueTADs > 5 )
  
  data.frame(
    obs_cat1=obs_cat1,
    obs_cat2=obs_cat2,
    obs_cat3=obs_cat3,
    obs_cat4=obs_cat4,
    rd_cat1=rd_cat1,
    rd_cat2=rd_cat2,
    rd_cat3=rd_cat3,
    rd_cat4=rd_cat4,
    stringsAsFactors = FALSE
  )
}))
rownames(all_cat_nbrUniqueTADs2) <- NULL


stopifnot(colSums(all_cat_nbrUniqueTADs) == colSums(all_cat_nbrUniqueTADs2))

## =>>> do stack barplot





# ratio rd/obs # unique TADs

dslist=all_ds_results[[1]]
famlist=dslist[[1]]

all_rd_over_obs_nbrUniqueTADs <- unlist(lapply(all_ds_results, function(dslist) lapply(dslist, function(famlist) {
  obs_uniqueTADs <- famlist[["nbrUniqueRegions"]]
  rd_uniqueTADs <- unlist(lapply(famlist[["random_results_data"]], function(x) x[["random_nbrUniqueRegions"]]))
  rd_uniqueTADs/obs_uniqueTADs
}))) 


plot(density(all_rd_over_obs_nbrUniqueTADs))
qt_05 <- quantile(all_rd_over_obs_nbrUniqueTADs, probs = 0.05)
qt_95 <- quantile(all_rd_over_obs_nbrUniqueTADs, probs = 0.95)
plot(density(all_rd_over_obs_nbrUniqueTADs[all_rd_over_obs_nbrUniqueTADs >= qt_05 & all_rd_over_obs_nbrUniqueTADs <= qt_95]))

s
mean(all_rd_over_obs_nbrUniqueTADs > 1)
mean(all_rd_over_obs_nbrUniqueTADs == 1)
mean(all_rd_over_obs_nbrUniqueTADs < 1)  # 0.01  #### faire un boxplot avec ces 3 catégories


# edges and unique TADs directly related:
#ratio #edges obs/@ edges rd
# # edges = # unqiue*(#unique-1) * 0.5














# all_obs_nbrCpts <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["nbrComponents"]])))
# 
# length(all_obs_nbrEdges)
# sum(all_obs_nbrEdges == 0)
# 
# length(all_obs_nbrCpts)
# sum(all_obs_nbrCpts == 0)
# 
# 

# 
# length(all_random_nbrEdges)
# sum(all_random_nbrEdges == 0)
# 
# 
# all_random_nbrCpts <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["random_results_dt"]]['random_nbrComponents',])))
# 
# length(all_random_nbrCpts)
# sum(all_random_nbrCpts == 0)
# 
# source("../Cancer_HiC_data_TAD_DA/utils_plot_fcts.R")
#     
# plot_multiDens(
#   list(
#     all_obs_nbrEdges=log10(0.01+all_obs_nbrEdges),
#     all_random_nbrEdges=log10(0.01+all_random_nbrEdges)
#   )
# )
# 
# plot_multiDens(
#   list(
#     all_obs_nbrCpts=log10(0.01+all_obs_nbrCpts),
#     all_random_nbrCpts=log10(0.01+all_random_nbrCpts)
#   )
# )
# 

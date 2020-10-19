
startTime <- Sys.time()

# Rscript sameFamInSameTAD_regionOcc_RANDOM.R

## !!! OUTPUT FOLDER:
# _RANDOM -> RANDOMMIDPOSSTRICT
# _RANDOMMIDPOSDISC -> RANDOMMIDPOSDISC
# _RANDOMMIDPOS -> RANDOMMIDPOS


# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity

set.seed(20200903) # NB forgot the first time I launched it

require(doMC)
require(foreach)
registerDoMC(40)
require(reshape2)
require(igraph)
require(ggpubr)
require(ggsci)

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

outFolder <- file.path("SAMEFAMINSAMETAD_REGIONOCC_RANDOMMIDPOS")
dir.create(outFolder, recursive = TRUE)

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]
#all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_hicds <- all_hicds[grepl("RANDOMMIDPOS_", all_hicds)]
# all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))


# all_ds <- unlist(sapply(all_hicds, function(x) file.path(x, list.files(file.path(pipFolder, x)))), use.names = FALSE)

hicds = "Barutcu_MCF-10A_40kb"
exprds = "TCGAbrca_lum_bas"
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]

buildData <- TRUE

# ds=all_ds[1]

# all_ds=all_ds[1]

if(buildData) {
  
  
  all_ds_results <- foreach(hicds = all_hicds) %dopar%{
    
    
    cat(paste0("... start: ", hicds,  "\n"))
    
    
    outFile <- file.path(outFolder, paste0(hicds, "_", "all_fam_results.Rdata"))
    if(file.exists(outFile)) {
      all_fam_results <- get(load(outFile))
      return(all_fam_results)
    }
    
    
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
      
      cat(paste0("... ", hicds, " - start fam: ", i_fam, "/", length(all_fams), "\n"))
      
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
    
# all_ds_results = get(load("SAMEFAMINSAMETAD_REGIONOCC/all_ds_results.Rdata"))

all_obs_ratioUniqueRegions <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["ratioUniqueRegions"]])))
stopifnot(all_obs_ratioUniqueRegions >= 0 & all_obs_ratioUniqueRegions <= 1)
all_random_ratioUniqueRegions <- unlist(lapply(all_ds_results, function(dslist) lapply(dslist, function(famlist) lapply(famlist[["random_results_data"]], function(x) x[["random_ratioUniqueRegions"]]))))
stopifnot(all_random_ratioUniqueRegions >= 0 & all_random_ratioUniqueRegions <= 1)

outFile <- file.path(outFolder,  paste0("allDS_ratioUniqueRegions_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(
    all_obs_ratioUniqueRegions=all_obs_ratioUniqueRegions,
    all_random_ratioUniqueRegions=all_random_ratioUniqueRegions
  ),
  plotTit ="ratioUniqueRegions"
)
mtext(side=3, text = paste0("all DS - n =", length(all_ds_results)), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))

# ratio of TADs with more than one gene
all_obs_ratioRegionsMultipleGenes <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) {
  mean(x[["regionOcc"]] > 1)
}))) 

all_random_ratioRegionsMultipleGenes <- unlist(lapply(all_ds_results, function(dslist) lapply(dslist, function(famlist) lapply(famlist[["random_results_data"]], function(x) {
  mean(x[["random_regionOcc"]] > 1)
}))))

outFile <- file.path(outFolder,  paste0("allDS_ratioRegionsMultipleGenes_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(
    all_obs_ratioRegionsMultipleGenes=all_obs_ratioRegionsMultipleGenes,
    all_random_ratioRegionsMultipleGenes=all_random_ratioRegionsMultipleGenes
  ),
  plotTit ="ratioRegionsMultipleGenes"
)
mtext(side=3, text = paste0("all DS - n =", length(all_ds_results)), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))


outFile <- file.path(outFolder,  paste0("allDS_ratioRegionsMultipleGenes_noZeros_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(
    all_obs_ratioRegionsMultipleGenes=all_obs_ratioRegionsMultipleGenes[all_obs_ratioRegionsMultipleGenes>0],
    all_random_ratioRegionsMultipleGenes=all_random_ratioRegionsMultipleGenes[all_random_ratioRegionsMultipleGenes>0]
  ),
  plotTit ="ratioRegionsMultipleGenes_noZeros"
)
mtext(side=3, text = paste0("all DS - n =", length(all_ds_results)), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))



#  barplot # with 1 unique Reg, 2-3 reg, 4-5, 5
all_obs_ratioUniqueRegions <- unlist(lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["ratioUniqueRegions"]])))

dslist=all_ds_results[[1]]
famlist=dslist[[1]]

# slower !
# all_cat_nbrUniqueTADs2 <- do.call(rbind, lapply(all_ds_results, function(dslist) do.call(rbind, lapply(dslist, function(famlist) {
#   obs_uniqueTADs <- famlist[["nbrUniqueRegions"]]
#   rd_uniqueTADs <- unlist(lapply(famlist[["random_results_data"]], function(x) x[["random_nbrUniqueRegions"]]))
#   
#   obs_cat1 <- sum(obs_uniqueTADs == 1)
#   obs_cat2 <- sum(obs_uniqueTADs == 2 | obs_uniqueTADs == 3 )
#   obs_cat3 <- sum(obs_uniqueTADs == 4 | obs_uniqueTADs == 5 )
#   obs_cat4 <- sum(obs_uniqueTADs > 5 )
#   
#   rd_cat1 <- sum(rd_uniqueTADs == 1)
#   rd_cat2 <- sum(rd_uniqueTADs == 2 | rd_uniqueTADs == 3 )
#   rd_cat3 <- sum(rd_uniqueTADs == 4 | rd_uniqueTADs == 5 )
#   rd_cat4 <- sum(rd_uniqueTADs > 5 )
#   
#   data.frame(
#     obs_cat1=obs_cat1,
#     obs_cat2=obs_cat2,
#     obs_cat3=obs_cat3,
#     obs_cat4=obs_cat4,
#     rd_cat1=rd_cat1,
#     rd_cat2=rd_cat2,
#     rd_cat3=rd_cat3,
#     rd_cat4=rd_cat4,
#     stringsAsFactors = FALSE
#   )
# }))) )
# rownames(all_cat_nbrUniqueTADs2) <- NULL
# stopifnot(colSums(all_cat_nbrUniqueTADs2) == colSums(all_cat_nbrUniqueTADs2))

# cat_labels <- c(
#   cat1 = "1 unique TAD",
#   cat2 = "2-3 unique TADs",
#   cat3 = "4-5 unique TADs",
#   cat4 = ">5 unique TADs"
#   )
cat_labels <- c(
  cat1 = "1",
  cat2 = "2-3",
  cat3 = "4-5",
  cat4 = ">5"
)



all_cat_nbrUniqueTADs <- do.call(rbind, lapply(all_ds_results, function(dslist) {
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
rownames(all_cat_nbrUniqueTADs) <- NULL

outFile <- file.path(outFolder, "all_cat_nbrUniqueTADs.Rdata")
save(all_cat_nbrUniqueTADs, file=outFile, version=2)

stopifnot(nrow(all_cat_nbrUniqueTADs) == length(all_ds_results))


## =>>> do stack barplot
plot_dt <- melt(all_cat_nbrUniqueTADs)
plot_dt$nbr_cat <- gsub(".+_(.+)", "\\1", plot_dt$variable)
plot_dt$data <- gsub("(.+)_.+", "\\1", plot_dt$variable)
plot_dt$data <- ifelse(plot_dt$data == "obs", "observed", 
                       ifelse(plot_dt$data == "rd", "random", NA))
plot_dt$nbr_cat_label <- cat_labels[plot_dt$nbr_cat]
stopifnot(!is.na(plot_dt))

tot_dt <- aggregate(value ~ data , FUN=sum, data=plot_dt)
tot_by_data <- setNames(tot_dt$value, tot_dt$data)
agg_dt <- aggregate(value ~ data + nbr_cat_label, FUN=sum, data=plot_dt)
agg_dt$value_ratio <- agg_dt$value/tot_by_data[agg_dt$data]

plotTit <- "# unique TADs by family (by occurence)"
subTit <- paste0("# DS = ", length(all_ds_results), "; # permut = ", nRandom)

agg_dt$value_log10 <- log10(agg_dt$value)

agg_dt$nbr_cat_label <- factor(agg_dt$nbr_cat_label, levels=cat_labels)
stopifnot(!is.na(agg_dt$nbr_cat_label))

p1_nbr <- ggbarplot(agg_dt, x="data", y="value_log10", fill="nbr_cat_label", 
                    xlab = "", ylab = "# unique TADs [log10]")+
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


plotTit <- "ratio unique TADs by family (by occurence)"

agg_dt$nbr_cat_label <- factor(agg_dt$nbr_cat_label, levels=cat_labels)
stopifnot(!is.na(agg_dt$nbr_cat_label))

p1_ratio <- ggbarplot(agg_dt, x="data", y="value_ratio", fill="nbr_cat_label",
                      xlab = "", ylab = "ratio unique TADs")+
  scale_fill_nejm()+
  labs(fill="") + 
  ggtitle(plotTit, subtitle=subTit)+
  theme(
    plot.title = element_text(size=16, face = "bold", hjust=0.5),
    plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
  )

outFile <- file.path(outFolder, paste0("nbrUniqueTADsByFam_ratio_by_cat_barplot.", plotType))
ggsave(p1_ratio, filename=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile,  "\n"))


# ratio rd/obs # unique TADs

dslist=all_ds_results[[1]]
famlist=dslist[[1]]

all_rd_over_obs_nbrUniqueTADs <- unlist(lapply(all_ds_results, function(dslist) lapply(dslist, function(famlist) {
  obs_uniqueTADs <- famlist[["nbrUniqueRegions"]]
  rd_uniqueTADs <- unlist(lapply(famlist[["random_results_data"]], function(x) x[["random_nbrUniqueRegions"]]))
  rd_uniqueTADs/obs_uniqueTADs
}))) 

outFile <- file.path(outFolder,  paste0("allDS_rd_over_obs_nbrUniqueTADs_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(all_rd_over_obs_nbrUniqueTADs), main="rd_over_obs_nbrUniqueTADs")
mtext(side=3, text = paste0("all DS - n =", length(all_ds_results)), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))


qt_05 <- quantile(all_rd_over_obs_nbrUniqueTADs, probs = 0.05)
qt_95 <- quantile(all_rd_over_obs_nbrUniqueTADs, probs = 0.95)
outFile <- file.path(outFolder,  paste0("allDS_rd_over_obs_nbrUniqueTADs_qtFilter_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(density(all_rd_over_obs_nbrUniqueTADs[all_rd_over_obs_nbrUniqueTADs >= qt_05 & all_rd_over_obs_nbrUniqueTADs <= qt_95]),
     main = "rd_over_obs_nbrUniqueTADs_qtFilter")
mtext(side=3, text = paste0("all DS - n =", length(all_ds_results)), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))

mean(all_rd_over_obs_nbrUniqueTADs > 1)
mean(all_rd_over_obs_nbrUniqueTADs == 1)
mean(all_rd_over_obs_nbrUniqueTADs < 1)  # 0.01  #### faire un boxplot avec ces 3 catégories


# edges and unique TADs directly related:
#ratio #edges obs/@ edges rd
# # edges = # unqiue*(#unique-1) * 0.5






cat(paste0(startTime, "\n", Sys.time(), "\n"))





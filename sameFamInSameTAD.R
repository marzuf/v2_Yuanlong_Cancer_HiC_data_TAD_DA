

# Rscript sameFamInSameTAD.R

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

aggFunc <- "mean"

plotType <- "svg"
myHeight <- 5
myWidth <- 7

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- file.path("SAMEFAMINSAMETAD_V2_ALLDS")
# outFolder <- file.path("SAMEFAMINSAMETAD")
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

buildData <- FALSE

ds=all_ds[1]

logOffset <- 0.01

# all_ds=all_ds[1]


### FOR CHECK !!!
setDir <- ""
hgnc_geneFamilyFile <- file.path(setDir, "/mnt/ed4/marie/family_data_2/hgnc_entrez_family.txt")
hgnc_geneFamilyDT <- read.delim(hgnc_geneFamilyFile, col.names=c("entrezID", "family"), header = F, stringsAsFactors = F)
hgnc_geneFamilyDT$entrezID <- as.character(hgnc_geneFamilyDT$entrezID)
hgnc_geneFamilyDT$family_short <- unlist(sapply(hgnc_geneFamilyDT$family, function(x) strsplit(x, "\\|")[[1]][1] ))
# any(duplicated(hgnc_geneFamilyDT$entrezID))


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
    tad_g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
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
    
    stopifnot(sameFam_dt$gene1 %in% g2t_dt$entrezID) 
    stopifnot(sameFam_dt$gene2 %in% g2t_dt$entrezID)
    stopifnot(sameFam_dt$gene1 %in% tad_g2t_dt$entrezID)
    stopifnot(sameFam_dt$gene2 %in% tad_g2t_dt$entrezID)
    
    
    sameFamSameTAD_dt <- merge(sameFam_dt, sameTAD_dt, all.x=TRUE, all.y=FALSE, by=c("gene1", "gene2"))
    
    all_genes <- familyDT$entrezID
    stopifnot(all_genes %in% names(entrezIDchromo))
    family_entrezIDchromo <- entrezIDchromo[names(entrezIDchromo) %in% all_genes]
    
    
    all_fams <- unique(sameFamSameTAD_dt$family)
    
    ##>> iterate here over families
    i_fam=1
    all_fam_results <- foreach(i_fam = 1:length(all_fams)) %dopar% {
      
      fam <- all_fams[i_fam]
      
      cat(paste0("... ", hicds, " - ", exprds, " - start fam: ", i_fam, "/", length(all_fams), "\n"))
      
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
    
# all_ds_results = get(load("SAMEFAMINSAMETAD/all_ds_results.Rdata"))

# lg1_1 = all_ds_results[["LG1_40kb/TCGAluad_mutKRAS_mutEGFR"]]
# lg1_2 = all_ds_results[["LG1_40kb/TCGAluad_norm_luad"]]
# lg1_1_genes <- lapply(lg1_1, function(x) x[["nbrGenes"]])
# lg1_2_genes <- lapply(lg1_2, function(x) x[["nbrGenes"]])
# all.equal(lg1_1_genes,lg1_2_genes)
# # TRUE
# all.equal(lg1_1[["nbrEdges"]], lg1_2[["nbrEdges"]])
# # TRUE
# all.equal(lg1_1[["nbrComponents"]], lg1_2[["nbrComponents"]])
# # TRUE

# singletons == TADs with regionOcc=1
# dt1 = all_ds_results[["Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas"]][["Tubulins"]]
# dt2 = all_ds_results2[["Barutcu_MCF-10A_40kb"]][["Tubulins"]]
# sum(dt2[["regionOcc"]] == 1) == dt1[["nbrSingletons"]]

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
          nbr_obs_dt$family <- gsub(".+/.+?\\.(.+)", "\\1", nbr_obs_dt$family_lab)  ### CORRECT HERE THERE DOTS IN THE FAMILY NAMES SHOULD USE INTERROGATION MARK 25.09.2020
          ### !!!! ADDED 25.09 IMPORTANT CHECK
          stopifnot(nbr_obs_dt$family %in% hgnc_geneFamilyDT$family_short)
          
          stopifnot(nbr_obs_dt$family != "")
          stopifnot(!is.na(nbr_obs_dt$family))
          
          nbr_rd_dt <- data.frame(family_lab=names(eval(parse(text=paste0("all_random_nbr", curr_var)))), 
                                       nbr=as.numeric(eval(parse(text=paste0("all_random_nbr", curr_var)))), 
                                       stringsAsFactors = FALSE)
          nbr_rd_dt$family <- gsub(".+/.+?\\.(.+)\\.result.+", "\\1", nbr_rd_dt$family_lab)  ### CORRECT HERE THERE DOTS IN THE FAMILY NAMES SHOULD USE INTERROGATION MARK 25.09.2020
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


}

# all_vars <- c("Edges", "Cpts", "Singletons")
all_vars <- c("Edges")
curr_var=all_vars[1]
for(curr_var in all_vars) {

  plotTit <- paste0("nbr", curr_var)
  subTit <- paste0("all DS - n =", length(all_ds_results))

  
  plot_list <- list(log10(logOffset+eval(parse(text = paste0("all_obs_nbr", curr_var)))),
                    log10(logOffset+eval(parse(text = paste0("all_random_nbr", curr_var)))))
  
  
  
  names(plot_list) <- c(paste0("all_obs_nbr",curr_var), paste0("all_random_nbr", curr_var) )

  outFile <- file.path(outFolder,  paste0("allDS_nbr", curr_var, "_density_", gsub("\\.", "", logOffset), ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(
      plot_list,
    plotTit = plotTit, 
    my_xlab = paste0("# ", curr_var, " [log10(+", logOffset, ")]")
  )
  mtext(side=3, text = subTit, font=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile,  "\n"))



  plotTit <- paste0("nbr", curr_var, "[log10(+", logOffset, ")]")
  subTit <- paste0("all DS - n =", length(all_ds_results))

  xlab <- "observed"
  ylab <- paste0("mean permut. (n=", nPermut, ")" )

  # outFile <- file.path(outFolder,  paste0("allDS_nbr", curr_var, "_meanPermut_vs_obs_densplot.", "svg"))
  outFile <- file.path(outFolder,  paste0("allDS_nbr", curr_var, "_meanPermut_vs_obs_densplot_", gsub("\\.", "", logOffset), ".", "svg"))
  do.call("svg", list(outFile, height=7, width=7))
  

  densplot(x=log10(eval(parse(text=paste0("all_obs_nbr", curr_var))) +logOffset),
           y=log10(eval(parse(text=paste0("all_random_meanNbr", curr_var))) +logOffset),
           xlab= xlab,ylab=ylab, main=plotTit
  )
  mtext(side=3, text = subTit, font=3)
  curve(1*x, lty=1, col="darkgrey", add = T)
  foo <- dev.off()
  cat(paste0("... written: ", outFile,  "\n"))


}


          
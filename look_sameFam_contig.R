

# Rscript look_sameFam_contig.R

# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity

# no need to have exprds

set.seed(20200903) # 

require(doMC)
require(foreach)
registerDoMC(40)
require(reshape2)
require(igraph)
require(ggpubr)
require(ggsci)

fontFamily <- "Hershey"

runFolder <- "."
pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")
familyVar <- "hgnc_family_short"

corMethod <- "pearson"
familyData <- "hgnc_family_short"

nRandom <- 5


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


outFolder <- file.path("LOOK_SAMEFAM_CONTIG")
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
    family_entrezIDchromo <- entrezIDchromo[names(entrezIDchromo) %in% all_genes] # used for the sampling
    
    all_fams <- unique(sameFamSameTAD_dt$family)

    
    ## !!! added here for gene rank
    hicds_gff_dt <- gff_dt[gff_dt$entrezID %in% names(family_entrezIDchromo),]
    hicds_gff_dt <- hicds_gff_dt[order(hicds_gff_dt$chromo, hicds_gff_dt$true_start, hicds_gff_dt$end),]

      
    # to avoid iterating over chromo -> when switching chromo, add a step of 2
    hicds_gff_dt$chromo <- as.character(hicds_gff_dt$chromo)
    
    # 17.10.2020: changed -> not needed, because i do this chromosome by chromosome
    # otherwise i could not select adjacent positions in the ranks
    hicds_gff_dt$gene_rank_offset <- sapply(hicds_gff_dt$chromo, function(x) which(unique(hicds_gff_dt$chromo) == x))
    hicds_gff_dt$gene_rank_offset <- hicds_gff_dt$gene_rank_offset - 1
    hicds_gff_dt$gene_rank_init <- 1:nrow(hicds_gff_dt)
    hicds_gff_dt$gene_rank <- hicds_gff_dt$gene_rank_offset + hicds_gff_dt$gene_rank_init
    # finally i dont need the break for the chromo because I sample by chromo
    # hicds_gff_dt[9121:9123,]
    # hicds_gff_dt[12078:12080,]
    stopifnot(sum(diff(hicds_gff_dt$gene_rank) != 1) == length(unique(hicds_gff_dt$chromo))-1)
    stopifnot(sum(diff(hicds_gff_dt$gene_rank) == 2) == length(unique(hicds_gff_dt$chromo))-1)
    # hicds_gff_dt$gene_rank <- 1:nrow(hicds_gff_dt)
    
    entrezID_geneRanks <- setNames(hicds_gff_dt$gene_rank, hicds_gff_dt$entrezID)
    
    stopifnot(setequal(names(family_entrezIDchromo), names(entrezID_geneRanks)))
    
    stopifnot(names(family_entrezIDchromo) %in% hicds_gff_dt$entrezID)
    
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
      
      ## should be done chromo by chromo !!!! -> ok for the look up
      true_ranks <- sort(entrezID_geneRanks[true_genes])
      stopifnot(length(true_ranks) == length(true_genes))
      names(true_ranks) <- NULL
      cont_gene_ranks <- diff(true_ranks)
      stopifnot(cont_gene_ranks > 0)
      rle_cont_genes <- rle(cont_gene_ranks)
      contig_tosample <- rle_cont_genes$lengths[rle_cont_genes$values==1]+1
      nbr_notcontig_tosample <- length(true_genes) - sum(contig_tosample)
      
      list(
        nbrGenes=nbrGenes,
        nbrSingletons=nbrSingletons,
        nbrEdges=nbrEdges,
        nbrComponents=nbrComponents,
        contig_tosample = contig_tosample
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
    
all_contig_to_sample <- lapply(all_ds_results, function(subl) lapply(subl, function(x) x[["contig_tosample"]]))

nDS <- length(all_contig_to_sample)

clusterSize <- unlist(all_contig_to_sample)
stopifnot(clusterSize > 1)
clusterSize_dt <- as.data.frame(clusterSize)

# for each dataset, how many contiguous clusters in each family
nClusters <- unlist(lapply(all_contig_to_sample, lengths))
nClusters_dt <- as.data.frame(nClusters)

# for each dataset, ratio of families that have at least 1 cluster
ratioFamWithClusters <- unlist(lapply(all_contig_to_sample, function(x) {
  nClust <- lengths(x)
  mean(nClust > 0)
  }))
ratioFamWithClusters_dt <- as.data.frame(ratioFamWithClusters)

all_vars <- c("clusterSize", "nClusters", "ratioFamWithClusters")
all_labs <- c("size (# genes) of the clusters", "# of clusters", "ratio families with at least 1 cluster")

plot_var <- "clusterSize"

subTit <- paste0("# DS = ", nDS)

i=1

for(i in 1:length(all_vars)) {
  
  plot_dt <- eval(parse(text=paste0(all_vars[i], "_dt")))
  
  plot_var <- all_vars[i]
  
  plotTit <- plot_var
  
  
  plot_dt[,paste0(plot_var, "_log10")] <- log10(plot_dt[,paste0(plot_var)])
  
  p3 <- ggdensity(plot_dt,
                  x =paste0(plot_var),
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0(all_labs[i]),
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  # color = "signif",
                  # fill = "signif",
                  palette = "jco"
  ) +
    ggtitle(plotTit, subtitle = subTit)+
    # scale_color_manual(values=my_cols)+
    # scale_fill_manual(values=my_cols)  +
    # labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
    # guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
      text = element_text(family=fontFamily),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
    )
  outFile <- file.path(outFolder, paste0(plot_var, "_all_ds_all_fams_density.", plotType))
  ggsave(p3, file=outFile, height=myHeight, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  if(plot_var=="clusterSize") {
    tmp_subTit <- paste0("# DS = ", nDS, " - clusterSize > 2 & < 20" )
    tmp_plot_dt <- plot_dt[plot_dt[,paste0(plot_var)] > 2 & plot_dt[,paste0(plot_var)] < 20,]
    p3 <- ggdensity(tmp_plot_dt,
                    x =paste0(plot_var),
                    y = "..density..",
                    # combine = TRUE,                  # Combine the 3 plots
                    xlab = paste0(all_labs[i]),
                    # add = "median",                  # Add median line.
                    rug = FALSE,                      # Add marginal rug
                    # color = "signif",
                    # fill = "signif",
                    palette = "jco"
    ) +
      ggtitle(plotTit, subtitle = tmp_subTit)+
      # scale_color_manual(values=my_cols)+
      # scale_fill_manual(values=my_cols)  +
      # labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
      # guides(color=FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      theme(
        text = element_text(family=fontFamily),
        panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x =  element_blank(),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
        axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
        axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
        plot.title = element_text(hjust=0.5, size = 16, face="bold"),
        plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
        legend.title = element_text(face="bold")
      )
    outFile <- file.path(outFolder, paste0(plot_var, "_all_ds_all_fams_density_subset.", plotType))
    ggsave(p3, file=outFile, height=myHeight, width=myWidth)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
  }
  
  
  if(plot_var=="nClusters") {
    tmp_subTit <- paste0("# DS = ", nDS, " - nClusters > 1 & < 20" )
    tmp_plot_dt <- plot_dt[plot_dt[,paste0(plot_var)] > 1 & plot_dt[,paste0(plot_var)] < 20,]
    p3 <- ggdensity(tmp_plot_dt,
                    x =paste0(plot_var),
                    y = "..density..",
                    # combine = TRUE,                  # Combine the 3 plots
                    xlab = paste0(all_labs[i]),
                    # add = "median",                  # Add median line.
                    rug = FALSE,                      # Add marginal rug
                    # color = "signif",
                    # fill = "signif",
                    palette = "jco"
    ) +
      ggtitle(plotTit, subtitle = tmp_subTit)+
      # scale_color_manual(values=my_cols)+
      # scale_fill_manual(values=my_cols)  +
      # labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
      # guides(color=FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      theme(
        text = element_text(family=fontFamily),
        panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x =  element_blank(),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
        axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
        axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
        plot.title = element_text(hjust=0.5, size = 16, face="bold"),
        plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
        legend.title = element_text(face="bold")
      )
    outFile <- file.path(outFolder, paste0(plot_var, "_all_ds_all_fams_density_subset.", plotType))
    ggsave(p3, file=outFile, height=myHeight, width=myWidth)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
  }
  
  
  
  
  
  p3 <- ggdensity(plot_dt,
                  x =paste0(plot_var, "_log10"),
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0(all_labs[i], " [log10]"),
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  # color = "signif",
                  # fill = "signif",
                  palette = "jco"
  ) +
    ggtitle(plotTit, subtitle = subTit)+
    # scale_color_manual(values=my_cols)+
    # scale_fill_manual(values=my_cols)  +
    # labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
    # guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
      text = element_text(family=fontFamily),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
    )
  outFile <- file.path(outFolder, paste0(plot_var, "_all_ds_all_fams_log10_density.", plotType))
  ggsave(p3, file=outFile, height=myHeight, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
  
  
  


}






startTime <- Sys.time()
cat(paste0("> Rscript tad_signif_nbr_enhancer.R\n"))

plotType <- "png"
myHeightGG  <- 7
myWidthGG <- 12


# Rscript tad_signif_nbr_enhancer.R

script_name <- "tad_signif_nbr_enhancer.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS=F
require(foreach)
require(doMC)
require(reshape2)
registerDoMC(ifelse(SSHFS, 2, 40))
require(ggplot2)


signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)
minOverlapBpRatio <- 0.8
minIntersectGenes <- 3


nTopConservedRegions <- 10

outFolder <- file.path("TAD_SIGNIF_NBR_ENHANCER")
dir.create(outFolder, recursive = TRUE)

buildTable <- FALSE


if(buildTable) {
  # enhancer data download from https://www.genecards.org/GeneHancer_version_4-4
  enhancerFile <- "GeneHancer_version_4-4.csv"
  stopifnot(file.exists(enhancerFile))
  enhancerDT <- read.delim(enhancerFile, stringsAsFactors = FALSE, sep="\t", header=TRUE)
  stopifnot(nrow(enhancerDT) > 0)
  
  
  final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
  stopifnot(file.exists(final_table_file))
  final_table_DT <- get(load(final_table_file))
  stopifnot(signif_column %in% colnames(final_table_DT))
  final_table_DT[,paste0(signifcol)] <- final_table_DT[,paste0(signif_column)] <= signifThresh
  
  final_table_DT$tad_id <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
  
  all_args <- c("", "norm_vs_tumor", "subtypes", "wt_vs_mut")
  args="subtypes"
  
  for(args in all_args) {
    if(args == "") {
      inFolder <- "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2"
      args <- ""
      file_prefix <- ""
    } else {
      inFolder <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", args[1])
      args <- args[1]
      file_prefix <- paste0(args, "_")
    } 
    stopifnot(dir.exists(inFolder))
    inFile <- file.path(inFolder, paste0(file_prefix,"conserved_regions_with_genes_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, "_df.txt"))
    stopifnot(file.exists(inFile))
    conserved_dt <- read.delim(inFile, header=TRUE, stringsAsFactors = FALSE)
    stopifnot(nrow(conserved_dt) >= nTopConservedRegions)
    conserved_dt <- conserved_dt[1:nTopConservedRegions,]
    stopifnot(nrow(conserved_dt) > 0)
    
    all_regions <- unlist(strsplit(conserved_dt$corresp_tads, split=","))
    new_col_name <- paste0(file_prefix, "top", nTopConservedRegions, "_conserved_signif")
    
    stopifnot(any(final_table_DT$tad_id %in% all_regions))
    
    final_table_DT[,paste0(new_col_name)] <- final_table_DT$tad_id %in% all_regions
  }
  
  nbrEnhancer_final_table_DT <- final_table_DT
  i=1
  nbrEnhancer_final_table_DT$nEnhancers <- foreach(i = 1:nrow(nbrEnhancer_final_table_DT), .combine='c') %dopar% {
    
    curr_region <- nbrEnhancer_final_table_DT$region[i]
    curr_chromo <- gsub("(chr.+)_.+", "\\1", curr_region)
    
    matching_enhancer_dt <- enhancerDT[enhancerDT$chrom == curr_chromo &
                                         enhancerDT$start >= nbrEnhancer_final_table_DT$start[i] &
                                         enhancerDT$start <= nbrEnhancer_final_table_DT$end[i],]
    
    nrow(matching_enhancer_dt)
    
    
  }
  
  outFile <- file.path(outFolder, "nbrEnhancer_final_table_DT.Rdata")
  save(nbrEnhancer_final_table_DT, file=outFile, version=2)
  
  cat(paste0("... written: ", outFile, "\n"))
  
}else { # end-if buildTable
  outFile="TAD_SIGNIF_NBR_ENHANCER/nbrEnhancer_final_table_DT.Rdata"
  outFile <- file.path(outFolder, "nbrEnhancer_final_table_DT.Rdata")
  nbrEnhancer_final_table_DT <- get(load(outFile))
  
}
keepCols <- 
keepCols <- c("tad_id", "nEnhancers", colnames(nbrEnhancer_final_table_DT)[grepl("conserved_signif", colnames(nbrEnhancer_final_table_DT))])

plot_DT <- nbrEnhancer_final_table_DT[,keepCols]
plot_DT_m <- melt(plot_DT, id=c("tad_id", "nEnhancers"))
plot_DT_m$plot_variable <- paste0(plot_DT_m$variable, plot_DT_m$value)


plot_DT_m$plot_variable <- gsub("conserved_signif", "conserved\nsignif", plot_DT_m$plot_variable)

p_var <-  ggplot(plot_DT_m, aes(x = plot_variable, y = nEnhancers)) + 
  geom_boxplot()+
  # facet_grid(~FDR_cut_off, switch="x") + 
  coord_cartesian(expand = FALSE) +
  # ggtitle(paste0(all_vars_tit[paste0(var_to_plot)], " by variable FDR cutoff"), subtitle = paste0("nDS = ", nDS))+
  ggtitle(paste0("# enhancers by TAD"), subtitle = paste0(""))+
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0("# enhancers"),
                     breaks = scales::pretty_breaks(n = 20))+
  # scale_fill_brewer(palette="YlOrRd")+
  # labs(fill  = "meanCorr type") +
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    strip.text.x = element_text(size = 10),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .3, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.text.x = element_text(angle=90, color="black", hjust=1,vjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )

outFile <- file.path(outFolder, paste0("nbr_enhancers_per_TAD_signif_notSignif_boxplot.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


##########################################################################################
##########################################################################################
##########################################################################################

cat("*** DONE - ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
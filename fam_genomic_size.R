

# Rscript fam_genomic_size.R


# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity


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

fontFamily <- "Hershey"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- file.path("FAM_GENOMIC_SIZE")
dir.create(outFolder, recursive = TRUE)

samFamFolder <- file.path(runFolder, "CREATE_SAME_FAMILY_SORTNODUP")
distFolder <- file.path(runFolder, "CREATE_DIST_SORTNODUP")

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
# all_hicds <- all_hicds[grepl("RANDOMMIDPOSSTRICT", all_hicds)]
# all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

# all_ds <- unlist(sapply(all_hicds, function(x) file.path(x, list.files(file.path(pipFolder, x)))), use.names = FALSE)

hicds = "Barutcu_MCF-10A_40kb"
# exprds = "TCGAbrca_lum_bas"
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]

buildData <- TRUE

# ds=all_ds[1]

# all_ds=all_ds[1]

if(buildData) {
  
  
  all_fam_dist_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    
    cat(paste0("... start: ", hicds,  "\n"))
    
    same_fam_dt <- get(load(file.path(samFamFolder, hicds, paste0(familyData, "_all_family_pairs.Rdata"))))
    # gene1  gene2   family
    # 1 10376  10381 Tubulins
    # 2 10376  10382 Tubulins
    # 3 10376  10383 Tubulins
    # 4 10376 112714 Tubulins
    same_fam_dt$gene1 <- as.character(same_fam_dt$gene1)
    same_fam_dt$gene2 <- as.character(same_fam_dt$gene2)
    stopifnot(same_fam_dt$gene1 < same_fam_dt$gene2)
    
    dist_dt <- get(load(file.path(distFolder, hicds, "all_dist_pairs.Rdata")))
    # gene1     gene2 chromo     dist
    # 1 100009667 100038246  chr10 33099818
    # 2 100009667     10006  chr10 42677402
    # 3 100009667 100118954  chr10 63950905
    # 4 100009667 100124332  chr10 55981525
    # 5 100009667 100125387  chr10 47856092
    dist_dt$gene1 <- as.character(dist_dt$gene1)
    dist_dt$gene2 <- as.character(dist_dt$gene2)
    stopifnot(dist_dt$gene1 < dist_dt$gene2)
    
    fam_dist_dt <- merge(same_fam_dt, dist_dt, by=c("gene1", "gene2"), all=FALSE)
    
    agg_mean_dt <- aggregate(dist ~ family, FUN=mean, data=fam_dist_dt)
    colnames(agg_mean_dt)[2] <- "mean_dist"
    agg_median_dt <- aggregate(dist ~ family, FUN=median, data=fam_dist_dt)
    colnames(agg_median_dt)[2] <- "median_dist"
    
    agg_dt <- merge(agg_mean_dt, agg_median_dt, by=c("family"), all=TRUE)
    stopifnot(!is.na(agg_dt))
    
    agg_dt$hicds <- hicds
    
    agg_dt
    
  }
  outFile <- file.path(outFolder, "all_fam_dist_dt.Rdata")
  save(all_fam_dist_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_fam_dist_dt.Rdata")
  # all_fam_dist_dt <- get(load("FAM_GENOMIC_SIZE/all_fam_dist_dt.Rdata"))
  all_fam_dist_dt <- get(load(outFile))
}

all_agg_funs <- c("mean", "median")
agg_fun="mean"

nDS <- length(unique(all_fam_dist_dt$hicds))

subTit <- paste0("# DS = ", nDS)

all_fam_dist_dt$log10_mean_dist <- log10(all_fam_dist_dt$mean_dist)
all_fam_dist_dt$log10_median_dist <- log10(all_fam_dist_dt$median_dist)

for(agg_fun in all_agg_funs) {
  
  plotTit <- paste0("dist. family genes (", agg_fun, ")")  
  
  
  p3 <- ggdensity(all_fam_dist_dt,
                  x =paste0("log10_", agg_fun, "_dist"),
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0(agg_fun, " dist. familiy gene pairs [bp log10]"),
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
  outFile <- file.path(outFolder, paste0("fam_gene_dist_", agg_fun, "_log10_density.", plotType))
  ggsave(p3, file=outFile, height=myHeight, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
}
    



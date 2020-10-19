#  Rscript cmp_samefaminsametad_v2_withRandom.R

outFolder <- "CMP_SAMEFAMINSAMETAD_V2_WITHRANDOM"
dir.create(outFolder, recursive = TRUE)

plotType <- "svg"
myHeight <- 6
myWidth <- 7
plotCex <- 1.2

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

require(reshape2)
logOffset <- 0.01


curr_var <- "Edges"
aggFunc <- "mean"
################################# observed
observed_data <- get(load("SAMEFAMINSAMETAD_V2_ALLDS/all_ds_results.Rdata"))
all_obs_nbrEdges <- unlist(lapply(observed_data, function(subl) lapply(subl, function(x) x[[paste0("nbr", curr_var)]])))

# aggregate by family
nbr_obs_dt <- data.frame(family_lab=names(eval(parse(text=paste0("all_obs_nbr", curr_var)))), 
                         nbr=as.numeric(eval(parse(text=paste0("all_obs_nbr", curr_var)))), stringsAsFactors = FALSE)
nbr_obs_dt$family <- gsub(".+/.+?\\.(.+)", "\\1", nbr_obs_dt$family_lab)  ### CORRECT HERE THERE DOTS IN THE FAMILY NAMES SHOULD USE INTERROGATION MARK 25.09.2020
nbr_obs_dt$ds <- file.path(gsub("(.+?)/.+", "\\1", nbr_obs_dt$family_lab), gsub(".+(TCGA.+?_.+?_.+?)\\..+", "\\1", nbr_obs_dt$family_lab))

stopifnot(nbr_obs_dt$family != "")
stopifnot(!is.na(nbr_obs_dt$family))

agg_obs_dt <- aggregate(nbr~family, data=nbr_obs_dt, FUN=aggFunc)
colnames(agg_obs_dt)[2] <- "nbr_obs"


################################# randommidposstrict
randomstrict_data <- get(load("SAMEFAMINSAMETAD_V2_ALLDS_RANDOM//all_ds_results.Rdata"))
all_randomstrict_nbrEdges <- unlist(lapply(randomstrict_data, function(subl) lapply(subl, function(x) x[[paste0("nbr", curr_var)]])))

# aggregate by family
nbr_randomstrict_dt <- data.frame(family_lab=names(eval(parse(text=paste0("all_randomstrict_nbr", curr_var)))), 
                         nbr=as.numeric(eval(parse(text=paste0("all_randomstrict_nbr", curr_var)))), stringsAsFactors = FALSE)
nbr_randomstrict_dt$family <- gsub(".+/.+?\\.(.+)", "\\1", nbr_randomstrict_dt$family_lab)  ### CORRECT HERE THERE DOTS IN THE FAMILY NAMES SHOULD USE INTERROGATION MARK 25.09.2020

nbr_randomstrict_dt$ds <- file.path(gsub("(.+?)/.+", "\\1", nbr_randomstrict_dt$family_lab), gsub(".+(TCGA.+?_.+?_.+?)\\..+", "\\1", nbr_randomstrict_dt$family_lab))

stopifnot( length(unique(nbr_obs_dt$ds)) == length(unique(nbr_randomstrict_dt$ds)))

stopifnot(nbr_randomstrict_dt$family != "")
stopifnot(!is.na(nbr_randomstrict_dt$family))

agg_randomstrict_dt <- aggregate(nbr~family, data=nbr_randomstrict_dt, FUN=aggFunc)
colnames(agg_randomstrict_dt)[2] <- "nbr_randommidposstrict"



################################# randommidposdisc

randomdisc_data <- get(load("SAMEFAMINSAMETAD_V2_ALLDS_RANDOM_RANDOMMIDPOSDISC///all_ds_results.Rdata"))
all_randomdisc_nbrEdges <- unlist(lapply(randomdisc_data, function(subl) lapply(subl, function(x) x[[paste0("nbr", curr_var)]])))

# aggregate by family
nbr_randomdisc_dt <- data.frame(family_lab=names(eval(parse(text=paste0("all_randomdisc_nbr", curr_var)))), 
                                  nbr=as.numeric(eval(parse(text=paste0("all_randomdisc_nbr", curr_var)))), stringsAsFactors = FALSE)
nbr_randomdisc_dt$family <- gsub(".+/.+?\\.(.+)", "\\1", nbr_randomdisc_dt$family_lab)  ### CORRECT HERE THERE DOTS IN THE FAMILY NAMES SHOULD USE INTERROGATION MARK 25.09.2020
nbr_randomdisc_dt$ds <- file.path(gsub("(.+?)/.+", "\\1", nbr_randomdisc_dt$family_lab), gsub(".+(TCGA.+?_.+?_.+?)\\..+", "\\1", nbr_randomdisc_dt$family_lab))

stopifnot( length(unique(nbr_obs_dt$ds)) == length(unique(nbr_randomdisc_dt$ds)))

stopifnot(nbr_randomdisc_dt$family != "")
stopifnot(!is.na(nbr_randomdisc_dt$family))

agg_randomdisc_dt <- aggregate(nbr~family, data=nbr_randomdisc_dt, FUN=aggFunc)
colnames(agg_randomdisc_dt)[2] <- "nbr_randommidposdisc"

##### PLOT

length(nbr_randomdisc_dt$ds)

all_agg_dt <- merge(agg_randomstrict_dt, merge(agg_obs_dt, agg_randomdisc_dt, by="family", all=TRUE), by="family", all=TRUE)


all_types <- c("randommidposstrict", "randommidposdisc")
rd_type <- c("randommidposdisc")

for(rd_type in all_types){
  
  plot_dt <- all_agg_dt[,c("family", "nbr_obs", paste0("nbr_", rd_type))]
  plot_dt <- na.omit(plot_dt)
  
  plotTit <- paste0("# of ", curr_var, " by family - obs. vs. ", toupper(rd_type))
  subTit <- paste0("aggreg. func = ", aggFunc, "; # families = ", length(unique(plot_dt$family)))
  
  plot_dt$nbr_obs_log10 <- log10(plot_dt$nbr_obs + logOffset)
  plot_dt[, paste0("nbr_", rd_type, "_log10")] <- log10(plot_dt[, paste0("nbr_", rd_type)] + logOffset)
  
  my_x <- plot_dt$nbr_obs_log10
  my_y <- plot_dt[, paste0("nbr_", rd_type, "_log10")]
  
  outFile <- file.path(outFolder,  paste0("allDS_nbr", curr_var, "_byFam_", rd_type, "_scatterplot_", gsub("\\.", "", logOffset), ".", plotType))
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



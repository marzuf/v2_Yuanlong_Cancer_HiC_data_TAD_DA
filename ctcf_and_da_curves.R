
# Rscript ctcf_and_da_curves.R

source("ctcf_da_utils.R")

library(ggplot2)

pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
all_obs_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_obs_exprds <- sapply(all_obs_hicds, function(x) list.files(file.path(pipFolder, x)))


# all_obs_exprds <- "TCGAluad_mutKRAS_mutEGFR"
# all_obs_hicds <- "ENCSR489OCU_NCI-H460_40kb"
# inFolder <- "CTCF_AND_DA_ALLDS_CHECK1"
inFolder <- "CTCF_AND_DA_ALLDS"


nBreaks <- 25
plotTypeGG <- "svg"
ggHeight <- 5
ggWidth <- 6

outFolder <- file.path("CTCF_AND_DA_CURVES")
dir.create(outFolder, recursive = TRUE)

inFile <- file.path(inFolder, "ctcf2tad_dt.Rdata")
ctcf2tad_dt <- get(load(inFile))
stopifnot(ctcf2tad_dt$hicds %in% all_obs_hicds)

inFile <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
final_dt <- get(load(inFile))
ds_final_dt <- final_dt
ds_final_dt <- final_dt[final_dt$hicds %in% unlist(all_obs_hicds) & final_dt$exprds %in% unlist(all_obs_exprds),  ]
colnames(ds_final_dt)[colnames(ds_final_dt) == "start"] <- "tad_start"
colnames(ds_final_dt)[colnames(ds_final_dt) == "end"] <- "tad_end"


both_dt <- merge(ds_final_dt, ctcf2tad_dt, by =c("hicds", "region"), all=FALSE)
stopifnot(nrow(both_dt) > 0)
both_dt$ctcf_midpos <- (both_dt$start + both_dt$end)/2
stopifnot(both_dt$end <= both_dt$tad_end)
stopifnot(both_dt$start >= both_dt$tad_start)

both_dt$relative_position <- (both_dt$ctcf_midpos - both_dt$tad_start)/(both_dt$tad_end - both_dt$tad_start)
stopifnot(both_dt$relative_position >= 0)
stopifnot(both_dt$relative_position <= 1)

rel_pos_breaks <- seq(from=0, to = 1, length.out=nBreaks)
rel_pos_labs <- get_fract_lab2(vect_values=both_dt$relative_position, range_levels = rel_pos_breaks)


rel_pos_levels <- gsub("^<=0$", "0",get_level_labs(rel_pos_breaks))

both_dt$rel_pos_lab <- rel_pos_labs
both_dt$rel_pos_lab <- factor(both_dt$rel_pos_lab, levels=rel_pos_levels)
stopifnot(!is.na(both_dt$rel_pos_lab))


head(both_dt)


############################################################################ by signif

all_to_plot <- list(
  c("relative_position", "length"), 
  c("MotifScore", "mean"), 
  c("MotifScore", "sum"), 
  c("ChipSeqScore", "mean"), 
  c("ChipSeqScore", "sum")
)

pthresh <- 0.01

both_dt$signif_lab <- ifelse(both_dt$adjPvalComb <= pthresh, "signif.", "not signif.")

my_theme <-theme_bw()+
  theme(
    axis.text.x = element_blank(),
    legend.position="top",
    # panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    # panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    # panel.background = element_rect(fill = "transparent"),
    # panel.grid.major.x =  element_blank(),
    # panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    plot.title = element_text(hjust=0.5, size = 14, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.text = element_text(size=12, hjust=0.5, vjust=0.5)
  )

  
for(i in 1:length(all_to_plot)) {
  
  
  curr_col <- all_to_plot[[i]][[1]]
  curr_fun <- all_to_plot[[i]][[2]]
  
  ########### by signif
  
  agg_byTAD_dt <- aggregate(as.formula(paste0(curr_col, " ~ hicds + exprds + region + rel_pos_lab + signif_lab")), data = both_dt, FUN=curr_fun)
  
  myTit <- paste0("# DS = ", length(unique(file.path(agg_byTAD_dt$hicds,agg_byTAD_dt$exprds))), 
                  "; # TADs = ", length(unique(file.path(agg_byTAD_dt$hicds,agg_byTAD_dt$exprds, agg_byTAD_dt$region ))))
  
  agg_all_dt <- aggregate(as.formula(paste0(curr_col, " ~ rel_pos_lab + signif_lab")), data=agg_byTAD_dt, FUN = "mean")
  
  p1 <- ggplot(data = agg_all_dt, aes_string(x="rel_pos_lab", y = paste0(curr_col), color="signif_lab", group="signif_lab")) + 
    labs(x="relative position in TAD", y=paste0(curr_fun, " agg. ", curr_col), color="")+
    ggtitle(myTit) +
    geom_line() + 
    scale_y_continuous(breaks=scales::pretty_breaks(n = 10))+
    my_theme
  
  outFile <- file.path(outFolder, paste0(curr_col, "_aggBy_", curr_fun, "_along_relPosInTAD_bySignif_lineplot.", plotTypeGG))
  ggsave(p1, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
  ########### by triplet class
  
  agg_byTAD_dt <- aggregate(as.formula(paste0(curr_col, " ~ hicds + exprds + region + rel_pos_lab + Triplet_class")), data = both_dt, FUN=curr_fun)
  
  myTit <- paste0("# DS = ", length(unique(file.path(agg_byTAD_dt$hicds,agg_byTAD_dt$exprds))), 
                  "; # TADs = ", length(unique(file.path(agg_byTAD_dt$hicds,agg_byTAD_dt$exprds, agg_byTAD_dt$region ))))
  
  agg_all_dt <- aggregate(as.formula(paste0(curr_col, " ~ rel_pos_lab + Triplet_class")), data=agg_byTAD_dt, FUN = "mean")
  
  p2 <- ggplot(data = agg_all_dt, aes_string(x="rel_pos_lab", y = paste0(curr_col), color="Triplet_class", group="Triplet_class")) + 
    labs(x="relative position in TAD", y=paste0(curr_fun, " agg. ", curr_col), color="")+
    ggtitle(myTit) +
    geom_line() + 
    scale_y_continuous(breaks=scales::pretty_breaks(n = 10))+
    my_theme
  
  
  outFile <- file.path(outFolder, paste0(curr_col, "_aggBy_", curr_fun, "_along_relPosInTAD_byTripletClass_lineplot.", plotTypeGG))
  ggsave(p2, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
}




# 
# 
# ############################################################################ by triplet class
# 
# curr_col <- "relative_position"
# curr_fun <- "length"
# 
# curr_col <- "MotifScore"
# curr_fun <- "mean"
# 
# pthresh <- 0.01
# 
# agg_byTAD_dt <- aggregate(as.formula(paste0(curr_col, " ~ hicds + exprds + region + rel_pos_lab + adjPvalComb")), data = both_dt, FUN=curr_fun)
# 
# agg_byTAD_dt$signif_lab <- ifelse(agg_byTAD_dt$adjPvalComb <= pthresh, "signif.", "not signif.")
# 
# agg_all_dt <- aggregate(as.formula(paste0(curr_col, " ~ rel_pos_lab + signif_lab")), data=agg_byTAD_dt, FUN = "mean")
# 
# ggplot(data = agg_all_dt, aes_string(x="rel_pos_lab", y = paste0(curr_col), color="signif_lab", group="signif_lab")) + 
#   geom_line() +
#   theme(axis.text.x = element_text(angle=90))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# plot(
#   x = as.numeric(agg_all_dt[,"rel_pos_lab"][agg_all_dt$signif_lab=="signif."]),
#   y = agg_all_dt[,paste0(curr_col)][agg_all_dt$signif_lab=="signif."], 
#   type="l"
# )
# 
# lines(
#   x = as.numeric(agg_all_dt[,"rel_pos_lab"][agg_all_dt$signif_lab=="not signif."]),
#   y = agg_all_dt[,paste0(curr_col)][agg_all_dt$signif_lab=="not signif."], 
#   type="l",
#   col="red"
# )
# 
# 
# 
# 
# 
# 
# 
# 
# # load("CTCF_AND_DA/ENCSR444WCZ_A549_40kb/TCGAluad_mutKRAS_mutEGFR/ctcf2tad_dt.Rdata")
# cat(paste0("# of CTCF BS:\t",nrow(ctcf2tad_dt), "\n"))
# cat(paste0("# of CTCF BS in TADs:\t",sum(!is.na(ctcf2tad_dt$region)), "\n"))
# cat(paste0("# of CTCF BS out of TADs:\t",sum(is.na(ctcf2tad_dt$region)), "\n"))
# 
# 
# tmpx <- unique(file.path(ds_final_dt$hicds))
# tmpy <- unique(file.path(ctcf2tad_dt$hicds))
# 
# cat(tmpy[!tmpy%in%tmpx], "\n")
# cat(tmpx[!tmpx%in%tmpy], "\n")
# 
# stopifnot(setequal(file.path(ds_final_dt$hicds),file.path(ctcf2tad_dt$hicds) ))
# 
# 
# merged_dt <- merge(ds_final_dt[,c("hicds", "exprds", "region", "meanLogFC", "meanCorr", "adjPvalComb", "tad_start", "tad_end")],
#                    ctcf2tad_dt,
#                    by=c("region", "hicds"), all=FALSE)
# tmp <- merged_dt$region
# tmp <- gsub("(.+)_.+", "\\1", tmp)
# stopifnot(tmp == merged_dt$chr)
# 
# merged_dt <- merged_dt[order(merged_dt$chr, merged_dt$start, merged_dt$end ),]
# # stopifnot(diff(merged_dt$start) >= 0) # not true because multiple chromo
# merged_dt$region <- as.character(merged_dt$region)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library("readxl")
# library(doMC)
# library(foreach)
# library(stringr)
# 
# 
# require(ggpubr)
# require(ggsci)
# 
# 
# pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
# 
# all_hicds <- list.files(pipFolder)
# all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
# all_obs_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
# all_obs_exprds <- sapply(all_obs_hicds, function(x) list.files(file.path(pipFolder, x)))
# 
# registerDoMC(40)
# # runFolder <- "../v2_Yuanlong_Cancer_HiC_data_TAD_DA_GM12878" #PIPELINE/OUTPUT_FOLDER/GM12878_40kb/TCGAluad_norm_luad/11sameNbr_runEmpPvalCombined/"
# # hicds <- "GM12878_40kb"
# 
# source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# plotType <- "png"
# myHeight <- myWidth <- 400
# plotCex <- 1.2
# 
# plotTypeGG <- "svg"
# ggHeight <- 6
# ggWidth <- 5
# 
# fontFamily <- "Hershey"
# 
# 
# runFolder <- "." 
# 
# 
# source("ctcf_da_utils.R")
# 
# final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
# ds_final_dt <- final_dt
# ds_final_dt <- final_dt[final_dt$hicds %in% unlist(all_obs_hicds) & final_dt$exprds %in% unlist(all_obs_exprds),  ]
# # ds_final_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds, ]
# # stopifnot(nrow(ds_final_dt) > 0)
# 
# buildTable <- TRUE
# 
# outFolder <- file.path("CTCF_AND_DA_ALLDS")
# dir.create(outFolder, recursive = TRUE)
# 
# init_ctcf_dt <- read_excel("13059_2020_2108_MOESM2_ESM.xlsx", sheet="CTCFs")
# init_ctcf_dt <- as.data.frame(init_ctcf_dt)
# init_ctcf_dt <- init_ctcf_dt[, 1:7]
# init_ctcf_dt$chr <- as.character(init_ctcf_dt$chr)
# 
# if(buildTable){
#   ctcf2tad_dt <- foreach(hicds = all_obs_hicds, .combine='rbind') %do% {
#     
#     cat(paste0("... start ", hicds, " \n"))
#     
#     tad_dt <- read.delim(file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt"), stringsAsFactors = FALSE, 
#                          header=F, col.names = c("chromo", "region", "start", "end"))
#     
#     ### KEEP ONLY TAD REGIONS
#     tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
#     stopifnot(!grepl("BOUND", tad_dt$region))
#     
#     # assign ctcf BS to tads
#     ctcf_dt <- init_ctcf_dt[init_ctcf_dt$chr %in% tad_dt$chromo,]
#     stopifnot(nrow(ctcf_dt) > 0)
#     stopifnot(is.numeric(ctcf_dt$start))
#     stopifnot(is.numeric(ctcf_dt$end))
#     stopifnot(ctcf_dt$start <= ctcf_dt$end)
#     ctcf_dt$region <- NA
#     ctcf_dt$hicds <- hicds
#     
#     i=1
#     i_ctcf2tad_dt <- foreach(i = 1:nrow(ctcf_dt), .combine='rbind') %dopar% {
#       
#       chr <- ctcf_dt$chr[i]
#       ctcf_start <- ctcf_dt$start[i]
#       ctcf_end <- ctcf_dt$end[i]
#       
#       subtad_dt <- tad_dt[tad_dt$chromo == chr,]
#       stopifnot(nrow(subtad_dt) > 0)
#       
#       # assign if start after tad start and end before tad end
#       test1 <- which(ctcf_start >= subtad_dt$start & ctcf_end <= subtad_dt$end)
#       test2 <- which(subtad_dt$start <= ctcf_start & subtad_dt$end >= ctcf_end)
#       stopifnot(test1==test2)
#       stopifnot(length(test1) == 0 | length(test1) == 1)
#       if(length(test1) == 1) {
#         ctcf_dt$region[i] <- subtad_dt$region[test1]
#       } else {
#         ctcf_dt$region[i] <- NA
#       }
#       ctcf_dt[i,]  
#     }
#     i_ctcf2tad_dt
#   } 
#   
#   outFile <- file.path(outFolder, "ctcf2tad_dt.Rdata")
#   save(ctcf2tad_dt, file=outFile, version=2)
#   cat(paste0("... written: ", outFile,"\n"))
#   
# }else {
# }

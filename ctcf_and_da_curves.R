
# Rscript ctcf_and_da_curves.R

get_fract_lab2 <- function(vect_values, range_levels) {
  stopifnot(is.numeric(vect_values))
  stopifnot(is.numeric(range_levels))
  range_levels <- sort(range_levels)
  my_cmd_first <- paste0("ifelse(vect_values <= ", range_levels[1], ",\"<=0\",")
  my_cmd_last <- paste0("ifelse(vect_values > ", range_levels[length(range_levels)], ",\">", round(range_levels[length(range_levels)],4),"\",NA")
  my_cmd_end <- paste0(rep(")", length(range_levels)+1), collapse="")
  my_cmd_mid <- paste0("ifelse(vect_values >",range_levels[1:(length(range_levels)-1)], " & vect_values <=", range_levels[2:(length(range_levels))],
                       ",\">",round(range_levels[1:(length(range_levels)-1)],4), " & <=", round(range_levels[2:(length(range_levels))],4),"\","
  )
  full_cmd <- paste0(my_cmd_first,
                     paste0(my_cmd_mid, collapse=""),
                     my_cmd_last,
                     my_cmd_end, collapse=",")
  return(eval(parse(text=full_cmd)))
}


inFolder <- "CTCF_AND_DA_ALLDS_CHECK1"

inFile <- file.path(inFolder, "ctcf2tad_dt.Rdata")
ctcf2tad_dt <- get(load(inFile))

inFile <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
ds_final_dt <- getload(inFile)
ds_final_dt <- final_dt[final_dt$hicds %in% unlist(all_obs_hicds) & final_dt$exprds %in% unlist(all_obs_exprds),  ]
colnames(ds_final_dt)[colnames(ds_final_dt) == "start"] <- "tad_start"
colnames(ds_final_dt)[colnames(ds_final_dt) == "end"] <- "tad_end"


both_dt <- merge(ds_final_dt, ctcf2tad_dt, by =c("hicds", "region"), all=FALSE)
both_dt$ctcf_midpos <- (both_dt$start + both_dt$end)/2
stopifnot(both_dt$end <= both_dt$tad_end)
stopifnot(both_dt$start >= both_dt$tad_start)

both_dt$relative_position <- (both_dt$ctcf_midpos - both_dt$tad_start)/(both_dt$tad_end - both_dt$tad_start)
stopifnot(both_dt$relative_position >= 0)
stopifnot(both_dt$relative_position <= 1)

rel_pos_breaks <- seq(from=0, to = 1, length.out=100)
rel_pos_labs <- get_fract_lab0(vect_values=both_dt$relative_position, range_levels = rel_pos_breaks)

rel_pos_breaks <- seq(from=0, to = 1, length.out=20)
rel_pos_labs <- get_fract_lab2(vect_values=both_dt$relative_position, range_levels = rel_pos_breaks)


rel_pos_levels <- gsub("^<=0$", "0",get_level_labs(rel_pos_breaks))

both_dt$rel_pos_lab <- rel_pos_labs
both_dt$rel_pos_lab <- factor(both_dt$rel_pos_lab, levels=rel_pos_levels)
stopifnot(!is.na(both_dt$rel_pos_lab))


head(both_dt)

curr_col <- "relative_position"
curr_fun <- "length"

pthresh <- 0.01

agg_byTAD_dt <- aggregate(as.formula(paste0(curr_col, " ~ hicds + exprds + region + rel_pos_lab + adjPvalComb")), data = both_dt, FUN=curr_fun)

agg_byTAD_dt$signif_lab <- ifelse(agg_byTAD_dt$adjPvalComb <= pthresh, "signif.", "not signif.")

agg_all_dt <- aggregate(as.formula(paste0(curr_col, " ~ rel_pos_lab + signif_lab")), data=agg_byTAD_dt, FUN = "mean")

ggplot(data = agg_all_dt, aes_string(x="rel_pos_lab", y = paste0(curr_col), color="signif_lab", group=1)) + 
  geom_line()



plot(
  x = as.numeric(agg_all_dt[,"rel_pos_lab"][agg_all_dt$signif_lab=="signif."]),
  y = agg_all_dt[,paste0(curr_col)][agg_all_dt$signif_lab=="signif."], 
  type="l"
)

lines(
  x = as.numeric(agg_all_dt[,"rel_pos_lab"][agg_all_dt$signif_lab=="not signif."]),
  y = agg_all_dt[,paste0(curr_col)][agg_all_dt$signif_lab=="not signif."], 
  type="l",
  col="red"
)








# load("CTCF_AND_DA/ENCSR444WCZ_A549_40kb/TCGAluad_mutKRAS_mutEGFR/ctcf2tad_dt.Rdata")
cat(paste0("# of CTCF BS:\t",nrow(ctcf2tad_dt), "\n"))
cat(paste0("# of CTCF BS in TADs:\t",sum(!is.na(ctcf2tad_dt$region)), "\n"))
cat(paste0("# of CTCF BS out of TADs:\t",sum(is.na(ctcf2tad_dt$region)), "\n"))


tmpx <- unique(file.path(ds_final_dt$hicds))
tmpy <- unique(file.path(ctcf2tad_dt$hicds))

cat(tmpy[!tmpy%in%tmpx], "\n")
cat(tmpx[!tmpx%in%tmpy], "\n")

stopifnot(setequal(file.path(ds_final_dt$hicds),file.path(ctcf2tad_dt$hicds) ))


merged_dt <- merge(ds_final_dt[,c("hicds", "exprds", "region", "meanLogFC", "meanCorr", "adjPvalComb", "tad_start", "tad_end")],
                   ctcf2tad_dt,
                   by=c("region", "hicds"), all=FALSE)
tmp <- merged_dt$region
tmp <- gsub("(.+)_.+", "\\1", tmp)
stopifnot(tmp == merged_dt$chr)

merged_dt <- merged_dt[order(merged_dt$chr, merged_dt$start, merged_dt$end ),]
# stopifnot(diff(merged_dt$start) >= 0) # not true because multiple chromo
merged_dt$region <- as.character(merged_dt$region)















library("readxl")
library(doMC)
library(foreach)
library(stringr)


require(ggpubr)
require(ggsci)


pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
all_obs_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_obs_exprds <- sapply(all_obs_hicds, function(x) list.files(file.path(pipFolder, x)))

registerDoMC(40)
# runFolder <- "../v2_Yuanlong_Cancer_HiC_data_TAD_DA_GM12878" #PIPELINE/OUTPUT_FOLDER/GM12878_40kb/TCGAluad_norm_luad/11sameNbr_runEmpPvalCombined/"
# hicds <- "GM12878_40kb"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
plotType <- "png"
myHeight <- myWidth <- 400
plotCex <- 1.2

plotTypeGG <- "svg"
ggHeight <- 6
ggWidth <- 5

fontFamily <- "Hershey"


runFolder <- "." 


source("ctcf_da_utils.R")

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
ds_final_dt <- final_dt
ds_final_dt <- final_dt[final_dt$hicds %in% unlist(all_obs_hicds) & final_dt$exprds %in% unlist(all_obs_exprds),  ]
# ds_final_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds, ]
# stopifnot(nrow(ds_final_dt) > 0)

buildTable <- TRUE

outFolder <- file.path("CTCF_AND_DA_ALLDS")
dir.create(outFolder, recursive = TRUE)

init_ctcf_dt <- read_excel("13059_2020_2108_MOESM2_ESM.xlsx", sheet="CTCFs")
init_ctcf_dt <- as.data.frame(init_ctcf_dt)
init_ctcf_dt <- init_ctcf_dt[, 1:7]
init_ctcf_dt$chr <- as.character(init_ctcf_dt$chr)

if(buildTable){
  ctcf2tad_dt <- foreach(hicds = all_obs_hicds, .combine='rbind') %do% {
    
    cat(paste0("... start ", hicds, " \n"))
    
    tad_dt <- read.delim(file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt"), stringsAsFactors = FALSE, 
                         header=F, col.names = c("chromo", "region", "start", "end"))
    
    ### KEEP ONLY TAD REGIONS
    tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
    stopifnot(!grepl("BOUND", tad_dt$region))
    
    # assign ctcf BS to tads
    ctcf_dt <- init_ctcf_dt[init_ctcf_dt$chr %in% tad_dt$chromo,]
    stopifnot(nrow(ctcf_dt) > 0)
    stopifnot(is.numeric(ctcf_dt$start))
    stopifnot(is.numeric(ctcf_dt$end))
    stopifnot(ctcf_dt$start <= ctcf_dt$end)
    ctcf_dt$region <- NA
    ctcf_dt$hicds <- hicds
    
    i=1
    i_ctcf2tad_dt <- foreach(i = 1:nrow(ctcf_dt), .combine='rbind') %dopar% {
      
      chr <- ctcf_dt$chr[i]
      ctcf_start <- ctcf_dt$start[i]
      ctcf_end <- ctcf_dt$end[i]
      
      subtad_dt <- tad_dt[tad_dt$chromo == chr,]
      stopifnot(nrow(subtad_dt) > 0)
      
      # assign if start after tad start and end before tad end
      test1 <- which(ctcf_start >= subtad_dt$start & ctcf_end <= subtad_dt$end)
      test2 <- which(subtad_dt$start <= ctcf_start & subtad_dt$end >= ctcf_end)
      stopifnot(test1==test2)
      stopifnot(length(test1) == 0 | length(test1) == 1)
      if(length(test1) == 1) {
        ctcf_dt$region[i] <- subtad_dt$region[test1]
      } else {
        ctcf_dt$region[i] <- NA
      }
      ctcf_dt[i,]  
    }
    i_ctcf2tad_dt
  } 
  
  outFile <- file.path(outFolder, "ctcf2tad_dt.Rdata")
  save(ctcf2tad_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile,"\n"))
  
}else {
}

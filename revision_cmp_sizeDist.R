require(foreach)

setDir <- "/media/electron"
setDir <- ""

hicds <- "LG1_40kb"
hierarchyFolder <- file.path(setDir, 
                             "/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/2.Results_review/2.Sup_tab_used_for_2nd_submission",
                             "Supplementary_Data_1_domain_hierarchies_127HiCmaps")

tad_files <- list.files(file.path(hicds, "FINAL_DOMAINS_WITH_SCORES"), pattern="final_domains_with_scores.txt", full.names = TRUE)
stopifnot(length(tad_files) == 22)

tad_file = tad_files[1]
all_tad_dt <- foreach(tad_file = tad_files, .combine='rbind') %dopar% { 
                      read.delim(tad_file,
                                 header=FALSE,col.names=c("chromo", "start", "end", "score"))
  }
all_tad_dt$size <- all_tad_dt$end - all_tad_dt$start + 1 
mean(all_tad_dt$size)

regions_dt <- read.delim("LG1_40kb/genes2tad/all_assigned_regions.txt", col.names=c("chromo", "region", "start", "end"))
regions_dt <- regions_dt[grepl("_TAD", regions_dt$region),]
regions_dt$size <- regions_dt$end - regions_dt$start + 1
mean(regions_dt$size)

stopifnot(mean(all_tad_dt$size) == mean(regions_dt$size))
stopifnot(length(setdiff(file.path(all_tad_dt$chromo, all_tad_dt$start, all_tad_dt$end),file.path(regions_dt$chromo, regions_dt$start, regions_dt$end))) == 0)

cptmt_file <- "Lung_tissue_1_binsize=40kb.bed"
cptmt_dt <- read.delim(file.path(hierarchyFolder, cptmt_file), header=F, skip=1, stringsAsFactors = FALSE,
                       col.names=c( "chr", "pos_start", "pos_end", "cptmt_label", "normalized_rank",".",
                                    "pos_start", "pos_end", "color", "1", "cptmt_score"))
stopifnot(cptmt_dt$pos_start.1 == cptmt_dt$pos_start)
stopifnot(cptmt_dt$pos_end.1 == cptmt_dt$pos_end)

cptmt_dt$size <- cptmt_dt$pos_end - cptmt_dt$pos_start + 1

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plot_multiDens(
  list(
    cptmt_size = log10(cptmt_dt$size),
    tad_size = log10(all_tad_dt$size)
  )
)

SSHFS <- FALSE
prefixDir <- ifelse(SSHFS, "/media/electron", "")
allData_file <- file.path(prefixDir, "/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/0.Scripts/Generate_TADs_for_Marie/bin_comp_uniq_hq_CELL_LINEs_list_debug.Rdata")
stopifnot(file.exists(allData_file))
allData_wrapped <- get(load(allData_file))
allData <- allData_wrapped$bin_comp_CELL_LINEs_list
data_info <- allData_wrapped$uniq_hq_CELL_LINE_info
all_ds <- names(allData)



folder_names <- setNames(all_ds, all_ds)
folder_names <- gsub("AWS_|Compendium_|mega_", "", folder_names)
stopifnot(!duplicated(folder_names))




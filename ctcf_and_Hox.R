
# Rscript ctcf_and_Hox.R


################# AJOUTER LE # TOT DANS LA LEGENDE !!!

# # source("ctcf_da_utils.R")
# 
# library(ggplot2)
# library(ggsci)

# pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
# 
# all_hicds <- list.files(pipFolder)
# all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
# all_obs_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
# all_obs_exprds <- sapply(all_obs_hicds, function(x) list.files(file.path(pipFolder, x)))


# all_obs_exprds <- "TCGAluad_mutKRAS_mutEGFR"
# all_obs_hicds <- "ENCSR489OCU_NCI-H460_40kb"
# inFolder <- "CTCF_AND_DA_ALLDS_CHECK1"
inFolder <- "CTCF_AND_DA_ALLDS"


# nBreaks <- 100
# step_breaks <- 1/nBreaks
# 
# plotTypeGG <- "svg"
# ggHeight <- 5
# ggWidth <- 6

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_mutKRAS_mutEGFR"

outFolder <- file.path("CTCF_AND_HOX")
dir.create(outFolder, recursive = TRUE)

inFile <- file.path(inFolder, "ctcf2tad_dt.Rdata")
ctcf2tad_dt <- get(load(inFile))
# stopifnot(ctcf2tad_dt$hicds %in% all_obs_hicds)

inFile <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
final_dt <- get(load(inFile))
ds_final_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds,]
stopifnot(nrow(ds_final_dt) > 0)
# ds_final_dt <- final_dt[final_dt$hicds %in% unlist(all_obs_hicds) & final_dt$exprds %in% unlist(all_obs_exprds),  ]
# colnames(ds_final_dt)[colnames(ds_final_dt) == "start"] <- "tad_start"
# colnames(ds_final_dt)[colnames(ds_final_dt) == "end"] <- "tad_end"

hox_final_dt <- ds_final_dt[grepl("^HOX", ds_final_dt$region_genes),]
# hicds                   exprds       region    start      end                        region_genes meanLogFC  meanCorr ratioDown meanLogFC_thresh_FDR0.1 meanLogFC_thresh_FDR0.2
# 76512  ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr17_TAD162 46520001 46720000 HOXB2,HOXB3,HOXB4,HOXB5,HOXB6,HOXB7 0.8883277 0.7062466         0                    1.25                     0.8
# 157612 ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR  chr7_TAD114 26920001 27200000 HOXA2,HOXA3,HOXA4,HOXA5,HOXA7,SKAP2 0.3390521 0.5667886         0                    1.25                     0.8

hox_tads <- as.character(hox_final_dt$region)
# add the hoxD tad
hox_tads <- c(hox_tads, "chr2_TAD643") # HOXD:chr2_TAD643 HOXA:chr7_TAD114 HOXB:chr17_TAD162

hox_ctcf_dt <- ctcf2tad_dt[ctcf2tad_dt$region %in% hox_tads &
              ctcf2tad_dt$hicds == hicds,
              ]

nMotifs <- aggregate(MotifScore ~ hicds + region, FUN=length, data=hox_ctcf_dt)
colnames(nMotifs)[colnames(nMotifs) == "MotifScore"] <- "nMotifs"
nMotifs

sumMotifScore <- aggregate(MotifScore ~ hicds + region, FUN=sum, data=hox_ctcf_dt)
colnames(sumMotifScore)[colnames(sumMotifScore) == "MotifScore"] <- "sumMotifScore"
sumMotifScore

meanMotifScore <- aggregate(MotifScore ~ hicds + region, FUN=mean, data=hox_ctcf_dt)
colnames(meanMotifScore)[colnames(meanMotifScore) == "MotifScore"] <- "meanMotifScore"
meanMotifScore

sumChipSeqScore <- aggregate(ChipSeqScore ~ hicds + region, FUN=sum, data=hox_ctcf_dt)
colnames(sumChipSeqScore)[colnames(sumChipSeqScore) == "ChipSeqScore"] <- "sumChipSeqScore"
sumChipSeqScore

meanChipSeqScore <- aggregate(ChipSeqScore ~ hicds + region, FUN=mean, data=hox_ctcf_dt)
colnames(meanChipSeqScore)[colnames(meanChipSeqScore) == "ChipSeqScore"] <- "meanChipSeqScore"
meanChipSeqScore

nMotifs_byOrientation <- aggregate(MotifScore ~ hicds + region+ orientation, FUN=length, data=hox_ctcf_dt)
colnames(nMotifs_byOrientation)[colnames(nMotifs_byOrientation) == "MotifScore"] <- "nMotifs"
nMotifs_byOrientation

nMotifs_byTC <- aggregate(MotifScore ~ hicds + region+ Triplet_class, FUN=length, data=hox_ctcf_dt)
colnames(nMotifs_byTC)[colnames(nMotifs_byTC) == "MotifScore"] <- "nMotifs"
nMotifs_byTC


inFile <- file.path(inFolder, "clustByTAD_dt.Rdata")
clust_dt <- get(load(inFile))
hox_clust_dt <- clust_dt[clust_dt$region %in% hox_tads &
                          clust_dt$hicds == hicds &
                          clust_dt$exprds == exprds,]

unique(hox_clust_dt[,c("hicds", "exprds", "region", "nConvergent")])

maxInClust <- aggregate(nInClust~hicds+exprds+region, data = hox_clust_dt, FUN=max)
colnames(maxInClust)[colnames(maxInClust) == "nInClust"] <- "maxInClust"
maxInClust

meanInClust <- aggregate(nInClust~hicds+exprds+region, data = hox_clust_dt, FUN=mean)
colnames(meanInClust)[colnames(meanInClust) == "nInClust"] <- "meanInClust"
meanInClust










# 
# 
# 
# both_dt <- merge(ds_final_dt, ctcf2tad_dt, by =c("hicds", "region"), all=FALSE)
# stopifnot(nrow(both_dt) > 0)
# both_dt$ctcf_midpos <- (both_dt$start + both_dt$end)/2
# stopifnot(both_dt$end <= both_dt$tad_end)
# stopifnot(both_dt$start >= both_dt$tad_start)
# 
# both_dt$relative_position <- (both_dt$ctcf_midpos - both_dt$tad_start)/(both_dt$tad_end - both_dt$tad_start)
# stopifnot(both_dt$relative_position >= 0)
# stopifnot(both_dt$relative_position <= 1)
# 
# # rel_pos_breaks <- seq(from=0, to = 1, length.out=nBreaks)
# # rel_pos_labs <- get_fract_lab0(vect_values=both_dt$relative_position, range_levels = rel_pos_breaks)
# # rel_pos_levels <- gsub("^<=0$", "0",get_level_labs(rel_pos_breaks))
# # both_dt$rel_pos_lab <- rel_pos_labs
# 
# both_dt$rel_pos_lab <- both_dt$relative_position %/% step_breaks
# both_dt$rel_pos_lab <- factor(both_dt$rel_pos_lab, levels=0:nBreaks)
# stopifnot(!is.na(both_dt$rel_pos_lab))
# 
# 
# head(both_dt)
# 
# 
# both_dt$rel_pos_lab2 <- both_dt$relative_position
# 
# 
# 
# 
# # marie@electron:/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA$ grep HOXD /mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt
# # 3239    chr2    176957532       176960666       GRCh37.p13      +       HOXD13
# # 3238    chr2    176964530       176965488       GRCh37.p13      +       HOXD12
# # 3237    chr2    176971721       176976563       GRCh37.p13      +       HOXD11
# # 3236    chr2    176981492       176984670       GRCh37.p13      +       HOXD10
# # 100506783       chr2    176984724       177001826       GRCh37.p13      -       HOXD-AS2
# # 3235    chr2    176987413       176989645       GRCh37.p13      +       HOXD9
# # 3234    chr2    176994422       176997423       GRCh37.p13      +       HOXD8
# # 3232    chr2    177001669       177043737       GRCh37.p13      +       HOXD3
# # 3233    chr2    177015122       177017951       GRCh37.p13      +       HOXD4
# # 401022  chr2    177037923       177053686       GRCh37.p13      -       HOXD-AS1
# # 3231    chr2    177053307       177055635       GRCh37.p13      +       HOXD1
# 
# 
# # marie@electron:/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA$ grep chr2_TAD643  ENCSR489OCU_NCI-H460_40kb_INITVERSION/genes2tad/all_genes_positions.txt 
# # 80856   chr2    176790410       176867018       chr2_TAD643
# # 344191  chr2    176944835       176948690       chr2_TAD643
# # 3239    chr2    176957532       176960666       chr2_TAD643
# # 3238    chr2    176964530       176965488       chr2_TAD643
# # 3237    chr2    176971721       176976563       chr2_TAD643
# # 3236    chr2    176981492       176984670       chr2_TAD643
# # 3235    chr2    176987413       176989645       chr2_TAD643
# # 100129455       chr2    176992529       176994359       chr2_TAD643
# # 3234    chr2    176994422       176997423       chr2_TAD643
# 
# # x = read.delim("ENCSR489OCU_NCI-H460_40kb/genes2tad/all_assigned_regions.txt", stringsAsFactors = FALSE, header=F)
# # x[x$V1 == "chr2" & x$V3 >= 176960000 & x$V4 <= 177060000,]
# # x2 = x[x$V1 == "chr2",]
# # 
# # y=read.delim("ENCSR489OCU_NCI-H460_40kb_CHECK141220/FINAL_DOMAINS/ENCSR489OCU_NCI-H460_chr2_YL_40kb_final_domains.txt", stringsAsFactors = FALSE, header=F)
# # y[y$V1 == "chr2" & y$V2 >= 176960000 & y$V3 <= 177060000,]
# # z=y[y$V1 == "chr2" & y$V2 >= 176960000 ,]
# # 
# # x2 = x[x$V1 == "chr2",]
# 
# 
# 

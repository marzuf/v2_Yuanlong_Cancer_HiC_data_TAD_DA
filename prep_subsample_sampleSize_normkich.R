# Rscript prep_subsample_sampleSize_normkich.R

set.seed(05042020)

dt1 <- get(load("../FIGURES_V5_YUANLONG/NSIGNIFGENES_NSIGNIFTADS_SUMSAMP/nSignif_dt_sumSample.Rdata"))
dt2 <- get(load("../FIGURES_V5_YUANLONG/NSIGNIFGENES_NSIGNIFTADS_MINSAMP/nSignif_dt_minSample.Rdata"))

dt <- merge(dt1, dt2, by=c("hicds", "exprds", "dataset"))
dt$otherSample <- dt$sumSample - dt$minNbrSample

hicds <- "ENCSR079VIJ_G401_40kb"
exprds <- "TCGAkich_norm_kich"

setDir <- "/media/electron"
setDir <- ""
sample1_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets/", exprds, "norm_ID.Rdata")
sample2_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets", exprds, "kich_ID.Rdata")

s1 <- get(load(sample1_file))
length(s1) # 25
s2 <- get(load(sample2_file))
length(s2) # 65
stopifnot(!any(s1 %in% s2))
stopifnot(!any(s2 %in% s1))

# 25 -> 10-> 15 -> 20
# 65 -> 26 -> 39 -> 52

all_sub_nbr1 <- c(10,15,20)
all_sub_nbr2 <- c(26, 39, 52)

stopifnot(all_sub_nbr1/all_sub_nbr2 == length(s1)/length(s2))

stopifnot(length(all_sub_nbr1) == length(all_sub_nbr2))

iSub = 1
for(iSub in 1:length(all_sub_nbr1)) {
  
  nSub2 <- all_sub_nbr2[iSub]
  
  sub_s2 <- sample(s2, size=nSub2, replace=F)
  stopifnot(!duplicated(sub_s2))
  stopifnot(length(sub_s2) == nSub2)
  
  new_file2 <- file.path(paste0(dirname(sample2_file), "_sub", iSub), basename(sample2_file))
  dir.create(dirname(new_file2))
  
  save(sub_s2, file = new_file2, version=2)
  cat(paste0("... written: ", new_file2, "\n"))
  
  ##########
  
  nSub1 <- all_sub_nbr1[iSub]
  
  sub_s1 <- sample(s1, size=nSub1, replace=F)
  stopifnot(!duplicated(sub_s1))
  stopifnot(length(sub_s1) == nSub1)
  
  new_file1 <- file.path(paste0(dirname(sample1_file), "_sub", iSub), basename(sample1_file))
  dir.create(dirname(new_file1))
  
  save(sub_s1, file = new_file1, version=2)
  cat(paste0("... written: ", new_file1, "\n"))
  
  
  stopifnot( length(sub_s1)/length(sub_s2) == length(s1)/length(s2))
  stopifnot( length(sub_s1)/length(sub_s2) == nSub1/nSub2)
  stopifnot( length(sub_s1) == nSub1)
  stopifnot( length(sub_s2) == nSub2)
  
}

# ... written: //mnt/ed4/marie/other_datasets/TCGAkich_norm_kich_sub1/kich_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets//TCGAkich_norm_kich_sub1/norm_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets/TCGAkich_norm_kich_sub2/kich_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets//TCGAkich_norm_kich_sub2/norm_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets/TCGAkich_norm_kich_sub3/kich_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets//TCGAkich_norm_kich_sub3/norm_ID.Rdata


# COPY TO EMULATE FOLDER 
# cp -r ENCSR079VIJ_G401_40kb ENCSR079VIJ_G401_RANDOMSUB1_40kb
# rm -rf ENCSR079VIJ_G401_RANDOMSUB1_40kb/FINAL_DOMAINS*
# cp -r ENCSR079VIJ_G401_40kb ENCSR079VIJ_G401_RANDOMSUB2_40kb
# rm -rf ENCSR079VIJ_G401_RANDOMSUB2_40kb/FINAL_DOMAINS*
# cp -r ENCSR079VIJ_G401_40kb ENCSR079VIJ_G401_RANDOMSUB3_40kb
# rm -rf ENCSR079VIJ_G401_RANDOMSUB3_40kb/FINAL_DOMAINS*


# then run_pipeline_sh but only step1 to create the setting file
# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB1_40kb TCGAkich_norm_kich
# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB2_40kb TCGAkich_norm_kich
# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB3_40kb TCGAkich_norm_kich

# change in the setting file:

# sample1_file <- "/mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf/lowInf_ID.Rdata"
# => sample1_file <- "/mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub20/lowInf_ID.Rdata"
# 
# sample2_file <- "/mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf/highInf_ID.Rdata"
# => sample2_file <- "/mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub20/highInf_ID.Rdata"

# then run_pipeline_sh but WITHOUT step1 to create the setting file
# all_sub_nbr1 <- c(10,15,20)
# all_sub_nbr2 <- c(26, 39, 52)

x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR079VIJ_G401_RANDOMSUB1_40kb/TCGAkich_norm_kich/0_prepGeneData/rna_fpkmDT.Rdata"))
stopifnot(ncol(x) == all_sub_nbr1[1] + all_sub_nbr2[1])

x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR079VIJ_G401_RANDOMSUB1_40kb/TCGAkich_norm_kich/1_runGeneDE/DE_rnaseqDT.Rdata"))
stopifnot(ncol(x) == all_sub_nbr1[1] + all_sub_nbr2[1])

x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR079VIJ_G401_RANDOMSUB2_40kb/TCGAkich_norm_kich/0_prepGeneData/rna_fpkmDT.Rdata"))
stopifnot(ncol(x) == all_sub_nbr1[2] + all_sub_nbr2[2])

x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR079VIJ_G401_RANDOMSUB2_40kb/TCGAkich_norm_kich/1_runGeneDE/DE_rnaseqDT.Rdata"))
stopifnot(ncol(x) == all_sub_nbr1[2] + all_sub_nbr2[2])


x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR079VIJ_G401_RANDOMSUB3_40kb/TCGAkich_norm_kich/0_prepGeneData/rna_fpkmDT.Rdata"))
stopifnot(ncol(x) == all_sub_nbr1[3] + all_sub_nbr2[3])


x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR079VIJ_G401_RANDOMSUB3_40kb/TCGAkich_norm_kich/1_runGeneDE/DE_rnaseqDT.Rdata"))
stopifnot(ncol(x) == all_sub_nbr1[3] + all_sub_nbr2[3])



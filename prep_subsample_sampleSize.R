# Rscript prep_subsample_sampleSize.R

set.seed(05042020)

dt1 <- get(load("../FIGURES_V5_YUANLONG/NSIGNIFGENES_NSIGNIFTADS_SUMSAMP/nSignif_dt_sumSample.Rdata"))
dt2 <- get(load("../FIGURES_V5_YUANLONG/NSIGNIFGENES_NSIGNIFTADS_MINSAMP/nSignif_dt_minSample.Rdata"))

dt <- merge(dt1, dt2, by=c("hicds", "exprds", "dataset"))
dt$otherSample <- dt$sumSample - dt$minNbrSample

setDir <- "/media/electron"
setDir <- ""
sample1_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf/lowInf_ID.Rdata")
sample2_file <- file.path(setDir, "/mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf/highInf_ID.Rdata")

s1 <- get(load(sample1_file))
s2 <- get(load(sample2_file))
stopifnot(!any(s1 %in% s2))
stopifnot(!any(s2 %in% s1))

all_sub_nbr <- seq(from=20, to = 100, by=20)

stopifnot(length(s1) == length(s2))

nSub=20
for(nSub in all_sub_nbr) {
  
  sub_s2 <- sample(s2, size=nSub, replace=F)
  stopifnot(!duplicated(sub_s2))
  stopifnot(length(sub_s2) == nSub)
  
  new_file2 <- file.path(paste0(dirname(sample2_file), "_sub", nSub), basename(sample2_file))
  dir.create(dirname(new_file2))
  
  save(sub_s2, file = new_file2, version=2)
  cat(paste0("... written: ", new_file2, "\n"))
  
  ##########
  
  sub_s1 <- sample(s1, size=nSub, replace=F)
  stopifnot(!duplicated(sub_s1))
  stopifnot(length(sub_s1) == nSub)
  
  new_file1 <- file.path(paste0(dirname(sample1_file), "_sub", nSub), basename(sample1_file))
  dir.create(dirname(new_file1))
  
  save(sub_s1, file = new_file1, version=2)
  cat(paste0("... written: ", new_file1, "\n"))
  
  
}

# ... written: //mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub20/highInf_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub20/lowInf_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub40/highInf_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub40/lowInf_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub60/highInf_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub60/lowInf_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub80/highInf_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub80/lowInf_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub100/highInf_ID.Rdata
# ... written: //mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub100/lowInf_ID.Rdata

# COPY TO EMULATE FOLDER 
# cp -r ENCSR312KHQ_SK-MEL-5_40kb ENCSR312KHQ_SK-MEL-5_RANDOMSUB40_40kb
# rm -rf ENCSR312KHQ_SK-MEL-5_RANDOMSUB40_40kb/FINAL_DOMAINS*

# then run_pipeline_sh but only step1 to create the setting file
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB40_40kb TCGAskcm_lowInf_highInf

# change in the setting file:

# sample1_file <- "/mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf/lowInf_ID.Rdata"
# => sample1_file <- "/mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub20/lowInf_ID.Rdata"
# 
# sample2_file <- "/mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf/highInf_ID.Rdata"
# => sample2_file <- "/mnt/ed4/marie/other_datasets/TCGAskcm_lowInf_highInf_sub20/highInf_ID.Rdata"

# then run_pipeline_sh but WITHOUT step1 to create the setting file

x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR312KHQ_SK-MEL-5_RANDOMSUB20_40kb/TCGAskcm_lowInf_highInf/0_prepGeneData/rna_fpkmDT.Rdata"))
stopifnot(ncol(x) == 2*20)

x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR312KHQ_SK-MEL-5_RANDOMSUB20_40kb/TCGAskcm_lowInf_highInf/1_runGeneDE/DE_rnaseqDT.Rdata"))
stopifnot(ncol(x) == 2*20)

x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR312KHQ_SK-MEL-5_RANDOMSUB40_40kb/TCGAskcm_lowInf_highInf/0_prepGeneData/rna_fpkmDT.Rdata"))
stopifnot(ncol(x) == 2*40)

x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR312KHQ_SK-MEL-5_RANDOMSUB40_40kb/TCGAskcm_lowInf_highInf/1_runGeneDE/DE_rnaseqDT.Rdata"))
stopifnot(ncol(x) == 2*40)


x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR312KHQ_SK-MEL-5_RANDOMSUB60_40kb/TCGAskcm_lowInf_highInf/0_prepGeneData/rna_fpkmDT.Rdata"))
stopifnot(ncol(x) == 2*60)


x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR312KHQ_SK-MEL-5_RANDOMSUB60_40kb/TCGAskcm_lowInf_highInf/1_runGeneDE/DE_rnaseqDT.Rdata"))
stopifnot(ncol(x) == 2*60)

x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR312KHQ_SK-MEL-5_RANDOMSUB80_40kb/TCGAskcm_lowInf_highInf/0_prepGeneData/rna_fpkmDT.Rdata"))
stopifnot(ncol(x) == 2*80)


x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR312KHQ_SK-MEL-5_RANDOMSUB80_40kb/TCGAskcm_lowInf_highInf/1_runGeneDE/DE_rnaseqDT.Rdata"))
stopifnot(ncol(x) == 2*80)

x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR312KHQ_SK-MEL-5_RANDOMSUB100_40kb/TCGAskcm_lowInf_highInf/0_prepGeneData/rna_fpkmDT.Rdata"))
stopifnot(ncol(x) == 2*100)

x = get(load("/media/electron//mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR312KHQ_SK-MEL-5_RANDOMSUB100_40kb/TCGAskcm_lowInf_highInf/1_runGeneDE/DE_rnaseqDT.Rdata"))
stopifnot(ncol(x) == 2*100)

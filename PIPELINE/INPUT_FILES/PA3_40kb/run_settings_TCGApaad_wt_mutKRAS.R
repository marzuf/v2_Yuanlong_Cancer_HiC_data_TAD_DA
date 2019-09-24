
    # > file written: Sat, 08 Dec 2018 12:56:24 +0100 

    # in this file, settings that are specific for a run on a dataset

    # gives path to output folder
    pipOutFold <- "OUTPUT_FOLDER/TCGApaad_wt_mutKRAS"

    # full path (starting with /mnt/...)
    # following format expected for the input
    # colnames = samplesID
    # rownames = geneID
    # !!! geneID are expected not difficulted

    # *************************************************************************************************************************
    # ************************************ SETTINGS FOR 0_prepGeneData
    # *************************************************************************************************************************

    # UPDATE 07.12.2018: for RSEM data, the "analog" FPKM file is provided separately (built in prepData)
    rna_fpkmDT_file <- "/mnt/ed4/marie/other_datasets/TCGApaad_wt_mutKRAS/fpkmDT.Rdata"

    rnaseqDT_file <- "/mnt/ed4/marie/other_datasets/TCGApaad_wt_mutKRAS/rnaseqDT_v2.Rdata"
    my_sep <- "\t"
    # input is Rdata or txt file ?
    # TRUE if the input is Rdata
    inRdata <- TRUE

    # can be ensemblID, entrezID, geneSymbol
    geneID_format <- "entrezID"
    stopifnot(geneID_format %in% c("ensemblID", "entrezID", "geneSymbol"))

    # are geneID rownames ? -> "rn" or numeric giving the column
    geneID_loc <- "rn"
    stopifnot(geneID_loc == "rn" | is.numeric(geneID_loc))

    removeDupGeneID <- TRUE

    # *************************************************************************************************************************
    # ************************************ SETTINGS FOR 1_runGeneDE
    # *************************************************************************************************************************

    # labels for conditions
    cond1 <- "wt"
    cond2 <- "mutKRAS"

    # path to sampleID for each condition - should be Rdata ( ! sample1 for cond1, sample2 for cond2 ! )
    sample1_file <- "/mnt/ed4/marie/other_datasets/TCGApaad_wt_mutKRAS/wt_ID.Rdata"
    sample2_file <- "/mnt/ed4/marie/other_datasets/TCGApaad_wt_mutKRAS/mutKRAS_ID.Rdata"

    minCpmRatio <- 20/888 

    inputDataType <- "RSEM"

    nCpu <- 20

    # number of permutations
nRandomPermut <- 10000
step8_for_permutGenes <- TRUE
step8_for_randomTADsFix <- FALSE
step8_for_randomTADsGaussian <- FALSE
step8_for_randomTADsShuffle <- FALSE
step14_for_randomTADsShuffle <- FALSE


# > file edited: Fri, 23 Aug 2019 09:45:54 +0200 

# path to output folder:
pipOutFold <- "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/PA3_40kb/TCGApaad_wt_mutKRAS"

# OVERWRITE THE DEFAULT SETTINGS FOR INPUT FILES - use TADs from the current Hi-C dataset 
TADpos_file <- paste0(setDir, "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PA3_40kb/genes2tad/all_assigned_regions.txt")
#chr1    chr1_TAD1       750001  1300000
#chr1    chr1_TAD2       2750001 3650000
#chr1    chr1_TAD3       3650001 4150000

gene2tadDT_file <- paste0(setDir, "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PA3_40kb/genes2tad/all_genes_positions.txt")
#LINC00115       chr1    761586  762902  chr1_TAD1
#FAM41C  chr1    803451  812283  chr1_TAD1
#SAMD11  chr1    860260  879955  chr1_TAD1
#NOC2L   chr1    879584  894689  chr1_TAD1

# overwrite main_settings.R: nCpu <- 25
nCpu <- 40

# *************************************************************************************************************************
# ************************************ SETTINGS FOR PERMUTATIONS (5#_, 8c_)
# *************************************************************************************************************************

# number of permutations
nRandomPermut <- 100000
gene2tadAssignMethod <- "maxOverlap"
nRandomPermutShuffle <- 100000
step8_for_permutGenes <- TRUE
step8_for_randomTADsFix <- FALSE
step8_for_randomTADsGaussian <- FALSE
step8_for_randomTADsShuffle <- FALSE
step14_for_randomTADsShuffle <- FALSE


            # added here 13.08.2019 to change the filtering of min. read counts
            rm(minCpmRatio)
            min_counts <- 5
            min_sampleRatio <- 0.8


            # to have compatible versions of Rdata
            options(save.defaults = list(version = 2))


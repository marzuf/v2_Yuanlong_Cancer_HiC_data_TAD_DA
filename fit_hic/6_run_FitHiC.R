

require(data.table)
require(FitHiC)

options(scipen=100)

script_name <- "6_run_FitHiC.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

# INPUT_AGG/KARPAS_DMSO/KARPAS_DMSO_chr1_noDS_merged_agg.txt 
# 0       1       1
# 0       14293   1
# 1       1       156
# 1       2       5


########################################################################################################
########################################################################################################
########################################################################################################
# FRAGSFILE This file stores the information about midpoints (or start indices) of the fragments. 
# 5 columns: 
#   1) chromosome name; 
#   2) will be ignored -> set to 0
#   3) midPoint
#   4) hitCount -> set to 1
#   5) will be ignored -> set to 0
# hitCount only used for filtering; set bins that need to be filtered to 0  [filtering can be combined with mappability threshold]
# FRAGSFILE Chromosome.Name 	Column.2 	Mid.Point 	Hit.Count 	Column.5
# chr1 	0 	20000 	1 	0
# chr1 	0 	60000 	1 	0


# INTERSFILE This file stores the information about interactions between fragment pairs. 
# 5 columns: 
  # 1) and 3) chromosome name
  # 2) and 4) midPoints of the fragment pair
  # 5) contact count between the two bins. 
# ! midpoints in this file match those in fragsfile
# ! use RAW contact counts and NOT the normalized counts 
# INTERSFILE Chromosome1.Name 	Mid.Point.1 	Chromosome2.Name 	Mid.Point.2 	Hit.Count
# chr1 	100020000 	chr1 	100100000 	201
# chr1 	100020000 	chr1 	100140000 	232

# Rscript 5_prep_FitHiC.R INPUT_AGG/KARPAS_DMSO KARPAS_DMSO noDS_merged_agg.txt chr1 10000 PREP_FITHIC


fithic_fragfile <- "PREP_FITHIC/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_downsample16_merged_agg_FitHiC_fragsfile.txt"
fithic_interfile <- file.path("PREP_FITHIC/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_downsample16_merged_agg_FitHiC_intersfile.txt")
ouptut_dir <- "FITHIC_OUTPUT"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
fithic_fragfile <- args[1]
fithic_interfile <- args[2]
output_dir <- args[3]

dir.create(output_dir, recursive=TRUE)

cat("> start FitHiC \n")
FitHiC(fragsfile=fithic_fragfile, intersfile=fithic_interfile, outdir=output_dir)
       
       
########################################################################################################
########################################################################################################
########################################################################################################
cat("*** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





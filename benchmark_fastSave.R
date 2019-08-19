library(microbenchmark)
# Rscript benchmark_fastSave.R

setDir=""

infile <- "PIPELINE/OUTPUT_FOLDER/Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/5_runPermutationsMedian/permutationsDT.Rdata"

inData <- get(load(infile))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

source(paste0(pipScriptDir, "/", "my_save_pigz.R")) # to use customed fastSave save.pigz()

stopifnot(file.exists(pigz_exec_path))


sink("benchmark_fast_save.txt")

cat(paste0(Sys.time(),"\n"))

microbenchmark(
  save = save(inData, file = "/tmp/benchmark1.Rdata"),
  saveRDS = saveRDS(inData, file = "/tmp/benchmark2.rds"),
  savePigz = my_save.pigz(inData, file = "/tmp/benchmark3.Rdata",  pigz_exec_path = pigz_exec_path),
  times = 5,
  control = list(warmup = 2)
) -> mb

mb

cat(paste0(Sys.time(),"\n"))

sink()

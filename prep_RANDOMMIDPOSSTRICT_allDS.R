
# Rscript prep_RANDOMMIDPOSSTRICT_allDS.R

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "check_diffTADs_randommidpos.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"

require(foreach)
require(doMC)
registerDoMC(40)

pipFolder<- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))


ratioSameTAD_thresh <- 0.8
discFolder <- file.path("CHECK_DIFFTADS_RANDOMMIDPOS_STRICT", ratioSameTAD_thresh)

all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[grepl("_RANDOMMIDPOS_", all_hicds) ]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))

hicds = all_hicds[1]
# all_hicds = all_hicds[1]
# all_hicds <-  all_hicds[grepl("460", all_hicds)]
all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  
  stopifnot(grepl("_RANDOMMIDPOS_", hicds))
  new_hicds <- gsub("_RANDOMMIDPOS_", "_RANDOMMIDPOSSTRICT_", hicds)
  
  mycmd <- paste("cp -r", hicds, new_hicds)
  cat(paste0("> ", mycmd, "\n"))
  system(mycmd)
  
  hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    cat("... start ", hicds, " - ", exprds, "\n")
    
    stopifnot(grepl("_RANDOMMIDPOS_", hicds))
    new_hicds <- gsub("_RANDOMMIDPOS_", "_RANDOMMIDPOSSTRICT_", hicds)
    
    newFolder <- file.path(pipOutFolder, new_hicds, exprds)
    dir.create(newFolder, recursive = TRUE)
    
    mycmd1 <- paste("cp -r", file.path(pipOutFolder, hicds, exprds, script0_name), newFolder)
    cat(paste0("> ", mycmd1, "\n"))
    system(mycmd1)
    
    cpfile1 <- paste(file.path(discFolder, hicds, exprds, "pipeline_regionList.Rdata"))
    cpfile2 <- paste(file.path(discFolder, hicds, exprds, "pipeline_geneList.Rdata"))
    stopifnot(file.exists(cpfile1))
    stopifnot(file.exists(cpfile2))
    
    
    todelfile <- file.path(newFolder, script0_name, "pipeline_regionList.Rdata")
    stopifnot(file.exists(todelfile))
    mycmd2 <- paste("rm -f ", todelfile)
    cat(paste0("> ", mycmd2, "\n"))
    system(mycmd2)
    mycmd2b <- paste("cp", cpfile1, todelfile)
    system(mycmd2b)
    
    todelfile2 <- file.path(newFolder, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(todelfile2))
    mycmd3 <- paste("rm -f ", todelfile2)
    cat(paste0("> ", mycmd3, "\n"))
    system(mycmd3)
    mycmd3b <- paste("cp", cpfile2, todelfile2)
    system(mycmd3b)
    
  }
}

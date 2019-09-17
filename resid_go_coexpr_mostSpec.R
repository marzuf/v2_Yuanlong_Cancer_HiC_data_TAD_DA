startTime <- Sys.time()
cat(paste0("> Rscript resid_go_coexpr_mostSpec.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

# Rscript resid_go_coexpr_mostSpec.R

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7 )
myWidth <- myHeight 

buildTable <- FALSE

corMethod <- "pearson"

### PREPARE GO DATA
ontologyType <- "BP"

## Bimap interface:
x <- org.Hs.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Convert to a list
all_genes_GO_list <- as.list(x[mapped_genes])

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

outFolder <- file.path("RESID_GO_COEXPR_MOST_SPEC")
dir.create(outFolder, recursive = TRUE)

mainFolder <- file.path(".")
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")

all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds


script0_name <- "0_prepGeneData"


####################################################################################################################################### >>> prepare the data
if(buildTable) {
  all_mostSpecGO_dist_coexpr_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    
    tadFile <- file.path("CREATE_SAME_TAD_SORTNODUP", hicds,  "all_TAD_pairs.Rdata")
    stopifnot(file.exists(tadFile))
    all_TAD_pairs <- get(load(tadFile))
    all_TAD_pairs$gene1 <- as.character(all_TAD_pairs$gene1)
    all_TAD_pairs$gene2 <- as.character(all_TAD_pairs$gene2)
    stopifnot(all_TAD_pairs$gene1 < all_TAD_pairs$gene2)
    
    all_TAD_pairs_init <- all_TAD_pairs
    rm(all_TAD_pairs)
    exprds = all_exprds[[paste0(hicds)]][1]
    
    distFile <- file.path("CREATE_DIST_SORTNODUP", hicds, "all_dist_pairs.Rdata")
    stopifnot(file.exists(distFile))
    all_dist_DT <- get(load(distFile))
    

    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start ", hicds, " - ", exprds, "\n")
      
      # > PREPARE geneList DATA
      dataset_pipDir <- file.path(pipFolder, hicds, exprds)
      stopifnot(dir.exists(dataset_pipDir))
      pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
      
      
      # > PREPARE sameTAD DATA
      # ! need to be done because TAD dt built at hicds level (not hicds/exprds level)
      all_TAD_pairs <- all_TAD_pairs_init[all_TAD_pairs_init$gene1 %in% pipeline_geneList & 
                                            all_TAD_pairs_init$gene2 %in% pipeline_geneList,]
      stopifnot(pipeline_geneList %in% all_TAD_pairs$gene1 | pipeline_geneList %in% all_TAD_pairs$gene2)      
      

      

      # > PREPARE coexpr DATA
      coexprFile <- file.path("CREATE_COEXPR_SORTNODUP",hicds, exprds, corMethod, "coexprDT.Rdata")
      stopifnot(file.exists(coexprFile))
      coexprDT <- get(load(coexprFile))
      coexprDT$gene1 <- as.character(coexprDT$gene1)
      coexprDT$gene2 <- as.character(coexprDT$gene2)
      stopifnot(pipeline_geneList %in% coexprDT$gene1 | pipeline_geneList %in% coexprDT$gene2)
      stopifnot(coexprDT$gene1 < coexprDT$gene2)
      coexprDT <- coexprDT[coexprDT$gene1 %in% pipeline_geneList & 
                             coexprDT$gene2 %in% pipeline_geneList,]
      
      # > PREPARE GO DATA
      goFile <- file.path("CREATE_SAME_GO_SORTNODUP", hicds, exprds, "all_GO_pairs.Rdata")
      stopifnot(file.exists(goFile))
      goDT <- get(load(goFile))
      
      
      merge_dt <- merge(goDT, coexprDT, all.x=TRUE, all.y=FALSE, by=c("gene1", "gene2"))
      stopifnot(!is.na(merge_dt$coexpr))
      stopifnot(nrow(merge_dt) == nrow(goDT))
      
      
      merge_dt2 <- merge(merge_dt,  all_TAD_pairs, all.x=TRUE, all.y=FALSE, by=c("gene1", "gene2"))
      stopifnot(!is.na(merge_dt2$coexpr))
      stopifnot(nrow(merge_dt2) == nrow(merge_dt2))
      
      dist_coexpr_tad_dt <- merge(merge_dt2, all_dist_DT, by=c("gene1", "gene2"), all.x=TRUE, all.y=FALSE)
      
      
      stopifnot(nrow(merge_dt2) == nrow(dist_coexpr_tad_dt))
      
      
      
        # for the genes located on different, distance not avaiable
        dist_coexpr_tad_dt <- dist_coexpr_tad_dt[!is.na(dist_coexpr_tad_dt$dist),]
        stopifnot(nrow(dist_coexpr_tad_dt) > 0)
        
        
        dist_coexpr_tad_dt$mapTAD <- ifelse(is.na(dist_coexpr_tad_dt$region), "diffTAD", "sameTAD")
        stopifnot(dist_coexpr_tad_dt$mapTAD[is.na(dist_coexpr_tad_dt$region)] == "diffTAD")
        
        dist_coexpr_tad_dt <- dist_coexpr_tad_dt[,c("gene1", "gene2", "GO", "coexpr", "dist", "mapTAD")]
        
        
        
        otherCols <- colnames(dist_coexpr_tad_dt)
        firstCols <- c("hicds", "exprds")
        dist_coexpr_tad_dt$hicds <- hicds
        dist_coexpr_tad_dt$exprds <- exprds
        dist_coexpr_tad_dt[,c(firstCols, otherCols)]

    } # end-foreach iterating over exprds
    exprds_dt
  }# end-foreach iterating over hicds
  
  outFile <- file.path(outFolder, paste0("all_mostSpecGO_dist_coexpr_dt.Rdata"))
  save(all_mostSpecGO_dist_coexpr_dt, file=outFile, version=2)  
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, paste0("all_mostSpecGO_dist_coexpr_dt.Rdata"))
  cat(paste0("... load data\t", Sys.time(), ""))
  all_mostSpecGO_dist_coexpr_dt <- get(load(outFile))
  cat(paste0("-", Sys.time(), "\n"))
}


outFile="RESID_GO_COEXPR_MOST_SPEC_oneDS/all_mostSpecGO_dist_coexpr_dt.Rdata"
all_mostSpecGO_dist_coexpr_dt=get(load(outFile))

# plot(coexpr~dist, data=all_mostSpecGO_dist_coexpr_dt)

outFile <- file.path(outFolder, paste0("diffTAD_sameTAD_sameGO_mostSpecGO_coexpr_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(coexpr~mapTAD, data=all_mostSpecGO_dist_coexpr_dt,
        ylab="coexpr.",
        xlab="",
        main=paste0("same GO coexpr. "))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



maxTAD_dist <- max(all_mostSpecGO_dist_coexpr_dt$dist[all_mostSpecGO_dist_coexpr_dt$mapTAD=="sameTAD"])
maxTAD_dist_dt <- all_mostSpecGO_dist_coexpr_dt[all_mostSpecGO_dist_coexpr_dt$dist <= maxTAD_dist,]

outFile <- file.path(outFolder, paste0("diffTAD_sameTAD_sameGO_mostSpecGO_coexpr_maxTADsizeLimited_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(coexpr~mapTAD, data=maxTAD_dist_dt,
        ylab="coexpr.",
        xlab="",
        main=paste0("same GO coexpr. "))
mtext(side=3, text=paste0("max TAD size limited"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


m1 <- lm(coexpr~dist, data=all_mostSpecGO_dist_coexpr_dt)
m1_resid <- residuals(m1)

all_mostSpecGO_dist_coexpr_dt$coexpr_dist_resid <- m1_resid

outFile <- file.path(outFolder, paste0("diffTAD_sameTAD_sameGO_mostSpecGO_coexprDistResid_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(coexpr_dist_resid~mapTAD, data=all_mostSpecGO_dist_coexpr_dt,
        ylab="coexpr~dist resid.",
        xlab="",
        main=paste0("same GO coexpr. ~ dist. resid "))
mtext(side=3, text=paste0("max TAD size limited"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))





######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


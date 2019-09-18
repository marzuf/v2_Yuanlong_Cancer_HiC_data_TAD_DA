script_name <- "go_specificity_geneLevel_tadLevel_mean.R"
options(scipen=100)

SSHFS=F

cat("> START ", script_name, "\n")

# Rscript go_specificity_geneLevel_tadLevel_mean.R 0.05 0.05
# Rscript go_specificity_geneLevel_tadLevel_mean.R 0.01 0.05

library(clusterProfiler)
library(ontologySimilarity)
library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)
library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)
library(ggpubr)

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"
# data(gene_GO_terms)
data(GO_IC)


setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))
startTime <- Sys.time()


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myHeightGG <- 7
myWidthGG <- myHeightGG*1.2


limmaCol <- "#00AFBB"
tadCol <- "#FC4E07"

buildTable <- TRUE

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

args <- commandArgs(trailingOnly = TRUE)

tads_signifThresh <- args[1]
genes_signifTresh <- args[2]

outFolder <- file.path("GO_SPECIFICITY_GENELEVEL_TADLEVEL_MEAN", paste0("tadPvalThresh", tads_signifThresh, "_genePvalThresh", genes_signifTresh))
dir.create(outFolder, recursive = TRUE)

logFile <- file.path(outFolder, "go_signif_geneLevel_tadLevel_logFile.txt")
if(buildTable) file.remove(logFile)


printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}

my_GO2TERM <- clusterProfiler:::get_GO2TERM_table()
my_GO2TERM$term <- paste0("GO_", toupper(gsub(" ", "_", my_GO2TERM$Term)))


printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}


go_signif_col <- "p.adjust"
go_signifThresh <- 0.05


txt <- paste0("... go_signif_col\t=\t", go_signif_col, "\n")
printAndLog(txt, logFile)
txt <- paste0("... go_signifThresh\t=\t", go_signifThresh, "\n")
printAndLog(txt, logFile)

file_suffix <- paste0("tadPvalThresh", tads_signifThresh, "_genePvalThresh", genes_signifTresh)

inFile <- file.path("GO_SIGNIF_GENELEVEL_TADLEVEL", paste0("tadPvalThresh", tads_signifThresh, "_genePvalThresh", genes_signifTresh), "all_go_enrich_list.Rdata")
stopifnot(file.exists(inFile))
all_go_enrich_list <- get(load(inFile))

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

if(buildTable) {
  
  hicds = all_hicds[1]
  all_go_ic_list <- foreach(hicds = all_hicds) %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_list <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      cat(paste0("... start: ", hicds, " - ", exprds, "\n"))
      
      go_signif_limma_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["limma_signif_enrich_resultDT"]]
      
      if(!is.null(go_signif_limma_dt)) {
        txt <- paste0(hicds, " - ", exprds, " - limma signif.: # annot. GO:\t", nrow(go_signif_limma_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_limma_dt <- go_signif_limma_dt[go_signif_limma_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - limma signif.: # signif. annot. GO:\t", nrow(go_signif_limma_dt), "\n")
        printAndLog(txt, logFile)
        
        nEnrichedGO_signif_limma <- nrow(go_signif_limma_dt)
        
        signif_limma_go_terms <- rownames(go_signif_limma_dt)
  
        signif_limma_go_ids <- my_GO2TERM$go_id[my_GO2TERM$term %in% signif_limma_go_terms]
        txt <- paste0(hicds, " - ", exprds, " - ... limma signif.: found GO ids matching GO terms:\t", length(signif_limma_go_ids), "/", length(signif_limma_go_terms), "\n")
        printAndLog(txt, logFile)
        
        signif_limma_go_ic <- GO_IC[names(GO_IC) %in% signif_limma_go_ids]
        txt <- paste0(hicds, " - ", exprds, " - ... limma signif.: found GO information content matching GO ids:\t", length(signif_limma_go_ic), "/", length(signif_limma_go_ids), "\n")
        printAndLog(txt, logFile)
        txt <- paste0(hicds, " - ", exprds, " - ... limma signif.: total retrieved data:\t", length(signif_limma_go_terms), " -> ", length(signif_limma_go_ic), "\n")
        printAndLog(txt, logFile)
        
        limma_signif_mean_go_ic <- mean(signif_limma_go_ic)
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... limma signif.: NULL \n")
        printAndLog(txt, logFile)
        
        limma_signif_mean_go_ic <- NULL
      }
      go_signif_tads_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["tad_signif_enrich_resultDT"]]
      
      if(!is.null(go_signif_tads_dt)) {
        txt <- paste0(hicds, " - ", exprds, " - TAD signif: # annot. GO:\t", nrow(go_signif_tads_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_tads_dt <- go_signif_tads_dt[go_signif_tads_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - TAD signif: # signif. annot. GO:\t", nrow(go_signif_tads_dt), "\n")
        printAndLog(txt, logFile)
        
        nEnrichedGO_signif_tads <- nrow(go_signif_tads_dt)
        
        signif_tad_go_terms <- rownames(go_signif_tads_dt)
        
        signif_tad_go_ids <- my_GO2TERM$go_id[my_GO2TERM$term %in% signif_tad_go_terms]
        txt <- paste0(hicds, " - ", exprds, " - ... TAD signif: found GO ids matching GO terms:\t", length(signif_tad_go_ids), "/", length(signif_tad_go_terms), "\n")
        printAndLog(txt, logFile)
        
        signif_tad_go_ic <- GO_IC[names(GO_IC) %in% signif_tad_go_ids]
        txt <- paste0(hicds, " - ", exprds, " - ... TAD signif: found GO information content matching GO ids:\t", length(signif_tad_go_ic), "/", length(signif_tad_go_ids), "\n")
        printAndLog(txt, logFile)
        txt <- paste0(hicds, " - ", exprds, " - ... TAD signif: total retrieved data:\t", length(signif_tad_go_terms), " -> ", length(signif_tad_go_ic), "\n")
        printAndLog(txt, logFile)
        
        tad_signif_mean_go_ic <- mean(signif_tad_go_ic)
        
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... TAD signif.: NULL \n")
        printAndLog(txt, logFile)
        tad_signif_mean_go_ic <- NULL
      }
      list(
        nEnrichedGO_signif_limma = nEnrichedGO_signif_limma,
        nEnrichedGO_signif_tads = nEnrichedGO_signif_tads,
        limma_signif_mean_go_ic = limma_signif_mean_go_ic,
        tad_signif_mean_go_ic = tad_signif_mean_go_ic
      )
    } # end-foreach iterating over exprds
    names(exprds_list) <- file.path(hicds, all_exprds[[paste0(hicds)]])
    exprds_list
  } # end-foreach iterating over hicds
  names(all_go_ic_list) <- all_hicds
  outFile <- file.path(outFolder, paste0("all_go_ic_list.Rdata"))
  save(all_go_ic_list, file = outFile, version=2)
  stopifnot(length(unlist(all_go_ic_list, recursive=FALSE)) == length(all_datasets))
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, paste0( "all_go_ic_list.Rdata"))
  cat("... load data\n")
  all_go_ic_list <- get(load(outFile))
}

all_nGO_signif_tads <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_tads"]])))
all_nGO_signif_limma <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_limma"]])))


all_ic_signif_tads <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["tad_signif_mean_go_ic"]])))
all_ic_signif_limma <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["limma_signif_mean_go_ic"]])))



######################################################################################## PLOT IC

plot_dt <- rbind(
  data.frame(
    signif_type="tad_signif",
    mean_ic = as.numeric(all_ic_signif_tads),
    stringsAsFactors = FALSE
  ),
  data.frame(
    signif_type="limma_signif",
    mean_ic = as.numeric(all_ic_signif_limma),
    stringsAsFactors = FALSE
  )
)

notNa_tads <- sum(!is.na(plot_dt$mean_ic[plot_dt$signif_type=="tad_signif"]))
notNa_limma <- sum(!is.na(plot_dt$mean_ic[plot_dt$signif_type=="limma_signif"]))

p <- ggdensity(plot_dt, 
               title = paste0("mean enriched GO IC"),
               subtitle=paste0("notNa_tads=", notNa_tads, "; notNa_limma=", notNa_limma),
               x = "mean_ic", 
               xlab="mean enriched GO IC",
               color = "signif_type", fill = "signif_type",
               add = "mean", rug = TRUE,
               palette = c(limmaCol, tadCol))

outFile <- file.path(outFolder, paste0("all_ds_enrichedGO_meanIC_", file_suffix, "_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))


p <- ggviolin(plot_dt, 
              title = paste0("mean enriched GO IC"),
              subtitle=paste0("notNa_tads=", notNa_tads, "; notNa_limma=", notNa_limma),
              x="signif_type",
               y = "mean_ic", 
               ylab="mean enriched GO IC",
              xlab="",
               color = "signif_type", fill = "signif_type",
               add = "mean", rug = TRUE,
               palette = c(limmaCol, tadCol))
p <- p + scale_y_continuous(breaks = seq(0, max(na.omit(plot_dt$mean_ic)), by = 20))

outFile <- file.path(outFolder, paste0("all_ds_enrichedGO_meanIC_", file_suffix, "_violin.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))



######################################################################################## PLOT IC


plot_dt <- rbind(
  data.frame(
    signif_type="tad_signif",
    nbr_GO = as.numeric(all_nGO_signif_tads),
    stringsAsFactors = FALSE
  ),
  data.frame(
    signif_type="limma_signif",
    nbr_GO = as.numeric(all_nGO_signif_limma),
    stringsAsFactors = FALSE
  )
)

notNa_tads <- sum(!is.na(plot_dt$nbr_GO[plot_dt$signif_type=="tad_signif"]))
notNa_limma <- sum(!is.na(plot_dt$nbr_GO[plot_dt$signif_type=="limma_signif"]))

p <- ggdensity(plot_dt, 
               title = paste0("# enriched GO"),
               subtitle=paste0("notNa_tads=", notNa_tads, "; notNa_limma=", notNa_limma),
               x = "nbr_GO", 
               xlab="# enriched GO",
               color = "signif_type", fill = "signif_type",
               add = "mean", rug = TRUE,
               palette = c(limmaCol, tadCol))


outFile <- file.path(outFolder, paste0("all_ds_nbr_enrichedGO_", file_suffix, "_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))


p <- ggviolin(plot_dt, 
              title = paste0("# enriched GO"),
              subtitle=paste0("notNa_tads=", notNa_tads, "; notNa_limma=", notNa_limma),
              x="signif_type",
              y = "nbr_GO", 
              ylab="# enriched GO",
              xlab="",
              color = "signif_type", fill = "signif_type",
              add = "mean", rug = TRUE,
              palette = c(limmaCol, tadCol))
p <- p + scale_y_continuous(breaks = seq(0, max(na.omit(plot_dt$nbr_GO)), by = 20))


outFile <- file.path(outFolder, paste0("all_ds_nbr_enrichedGO_", file_suffix, "_violin.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




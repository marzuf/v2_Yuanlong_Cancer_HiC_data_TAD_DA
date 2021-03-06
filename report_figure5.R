script_name <- "report_figure5.R"
options(scipen=100)

SSHFS=F

cat("> START ", script_name, "\n")

# Rscript report_figure5.R 0.01 0.05

# Rscript report_figure5.R 0.05 0.05
# Rscript report_figure5.R 0.01 0.01

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


col1 <- get_palette("Dark2", 5)[1]
col2 <- get_palette("Dark2", 5)[2]
col3 <- get_palette("Dark2", 5)[3]
col4 <- get_palette("Dark2", 5)[4]
col5 <- get_palette("Dark2", 5)[5]


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
genes_signifThresh <- args[2]

outFolder <- file.path("REPORT_FIGURE5", paste0("tadPvalThresh", tads_signifThresh, "_genePvalThresh", genes_signifThresh))
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

file_suffix <- paste0("tadPvalThresh", tads_signifThresh, "_genePvalThresh", genes_signifThresh)

inFile <- file.path("GO_SIGNIF_GENELEVEL_TADLEVEL_INTERSECTDIFF", paste0("tadPvalThresh", tads_signifThresh, "_genePvalThresh", genes_signifThresh), "all_go_enrich_list.Rdata")
stopifnot(file.exists(inFile))
all_go_enrich_list <- get(load(inFile))

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

if(buildTable) {
  
  hicds = all_hicds[1]
  all_go_ic_list <- foreach(hicds = all_hicds) %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_list <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      cat(paste0("... start: ", hicds, " - ", exprds, "\n"))
      
      
      ############# intersect
      
      go_signif_intersect_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["intersect_signif_enrich_resultDT"]]
      
      if(!is.null(go_signif_intersect_dt)) {
        txt <- paste0(hicds, " - ", exprds, " - intersect signif.: # annot. GO:\t", nrow(go_signif_intersect_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_intersect_dt <- go_signif_intersect_dt[go_signif_intersect_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - intersect signif.: # signif. annot. GO:\t", nrow(go_signif_intersect_dt), "\n")
        printAndLog(txt, logFile)
        
        nEnrichedGO_signif_intersect <- nrow(go_signif_intersect_dt)
        
        signif_intersect_go_terms <- rownames(go_signif_intersect_dt)
        
        signif_intersect_go_ids <- my_GO2TERM$go_id[my_GO2TERM$term %in% signif_intersect_go_terms]
        txt <- paste0(hicds, " - ", exprds, " - ... intersect signif.: found GO ids matching GO terms:\t", length(signif_intersect_go_ids), "/", length(signif_intersect_go_terms), "\n")
        printAndLog(txt, logFile)
        
        signif_intersect_go_ic <- GO_IC[names(GO_IC) %in% signif_intersect_go_ids]
        txt <- paste0(hicds, " - ", exprds, " - ... intersect signif.: found GO information content matching GO ids:\t", length(signif_intersect_go_ic), "/", length(signif_intersect_go_ids), "\n")
        printAndLog(txt, logFile)
        txt <- paste0(hicds, " - ", exprds, " - ... intersect signif.: total retrieved data:\t", length(signif_intersect_go_terms), " -> ", length(signif_intersect_go_ic), "\n")
        printAndLog(txt, logFile)
        
        
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... intersect signif.: NULL \n")
        printAndLog(txt, logFile)
        
        nEnrichedGO_signif_intersect <- NA #changed 12.10
        
        signif_intersect_go_ic <- NA #changed 12.10
      }
      
      
      ############# limma only
      
      go_signif_limmaOnly_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["limmaOnly_signif_enrich_resultDT"]]
      
      if(!is.null(go_signif_limmaOnly_dt)) {
        txt <- paste0(hicds, " - ", exprds, " - limma only signif.: # annot. GO:\t", nrow(go_signif_limmaOnly_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_limmaOnly_dt <- go_signif_limmaOnly_dt[go_signif_limmaOnly_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - limma only signif.: # signif. annot. GO:\t", nrow(go_signif_limmaOnly_dt), "\n")
        printAndLog(txt, logFile)
        
        nEnrichedGO_signif_limmaOnly <- nrow(go_signif_limmaOnly_dt)
        
        signif_limmaOnly_go_terms <- rownames(go_signif_limmaOnly_dt)
        
        signif_limmaOnly_go_ids <- my_GO2TERM$go_id[my_GO2TERM$term %in% signif_limmaOnly_go_terms]
        txt <- paste0(hicds, " - ", exprds, " - ... limma only signif.: found GO ids matching GO terms:\t", length(signif_limmaOnly_go_ids), "/", length(signif_limmaOnly_go_terms), "\n")
        printAndLog(txt, logFile)
        
        signif_limmaOnly_go_ic <- GO_IC[names(GO_IC) %in% signif_limmaOnly_go_ids]
        txt <- paste0(hicds, " - ", exprds, " - ... limma only signif.: found GO information content matching GO ids:\t", length(signif_limmaOnly_go_ic), "/", length(signif_limmaOnly_go_ids), "\n")
        printAndLog(txt, logFile)
        txt <- paste0(hicds, " - ", exprds, " - ... limma only signif.: total retrieved data:\t", length(signif_limmaOnly_go_terms), " -> ", length(signif_limmaOnly_go_ic), "\n")
        printAndLog(txt, logFile)
        
        
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... limma only signif.: NULL \n")
        printAndLog(txt, logFile)
        
        signif_limmaOnly_go_ic <- NA #changed 12.10
        
        nEnrichedGO_signif_limmaOnly <- NA #changed 12.10
      }
      
      
      
      
      ############# TAD only
      
      
      
      go_signif_tadsOnly_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["tadsOnly_signif_enrich_resultDT"]]
      
      if(!is.null(go_signif_tadsOnly_dt)) {
        txt <- paste0(hicds, " - ", exprds, " - TAD only signif: # annot. GO:\t", nrow(go_signif_tadsOnly_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_tadsOnly_dt <- go_signif_tadsOnly_dt[go_signif_tadsOnly_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - TAD only signif: # signif. annot. GO:\t", nrow(go_signif_tadsOnly_dt), "\n")
        printAndLog(txt, logFile)
        
        nEnrichedGO_signif_tadsOnly <- nrow(go_signif_tadsOnly_dt)
        
        signif_tadOnly_go_terms <- rownames(go_signif_tadsOnly_dt)
        
        signif_tadOnly_go_ids <- my_GO2TERM$go_id[my_GO2TERM$term %in% signif_tadOnly_go_terms]
        txt <- paste0(hicds, " - ", exprds, " - ... TAD only signif: found GO ids matching GO terms:\t", length(signif_tadOnly_go_ids), "/", length(signif_tadOnly_go_terms), "\n")
        printAndLog(txt, logFile)
        
        signif_tadOnly_go_ic <- GO_IC[names(GO_IC) %in% signif_tadOnly_go_ids]
        txt <- paste0(hicds, " - ", exprds, " - ... TAD only signif: found GO information content matching GO ids:\t", length(signif_tadOnly_go_ic), "/", length(signif_tadOnly_go_ids), "\n")
        printAndLog(txt, logFile)
        txt <- paste0(hicds, " - ", exprds, " - ... TAD only signif: total retrieved data:\t", length(signif_tadOnly_go_terms), " -> ", length(signif_tadOnly_go_ic), "\n")
        printAndLog(txt, logFile)
        
        
        
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... TAD only signif.: NULL \n")
        printAndLog(txt, logFile)
        signif_tadOnly_go_ic <- NA #changed 12.10
        
        nEnrichedGO_signif_tadsOnly <- NA #changed 12.10
        
      }
      
      
      
      ############# limma
      
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
        
        
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... limma signif.: NULL \n")
        printAndLog(txt, logFile)
        
        signif_limma_go_ic <- NA #changed 12.10
        
        nEnrichedGO_signif_limma <- NA #changed 12.10
        
        
      }
      
      ############# TAD
      
      
      
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
        
        
        
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... TAD signif.: NULL \n")
        printAndLog(txt, logFile)
        signif_tad_go_ic <- NA #changed 12.10
        
        nEnrichedGO_signif_tads <- NA #changed 12.10
        
      }
      
      
      
      list(
        nEnrichedGO_signif_limma = nEnrichedGO_signif_limma,
        nEnrichedGO_signif_limmaOnly = nEnrichedGO_signif_limmaOnly,
        nEnrichedGO_signif_tads = nEnrichedGO_signif_tads,
        nEnrichedGO_signif_tadsOnly = nEnrichedGO_signif_tadsOnly,
        nEnrichedGO_signif_intersect = nEnrichedGO_signif_intersect,
        signif_limma_go_ic = signif_limma_go_ic,
        signif_limmaOnly_go_ic = signif_limmaOnly_go_ic,
        signif_tad_go_ic = signif_tad_go_ic,
        signif_tadOnly_go_ic = signif_tadOnly_go_ic,
        signif_intersect_go_ic = signif_intersect_go_ic
      )
      
      
    } # end-foreach iterating over exprds
    names(exprds_list) <- file.path(hicds, all_exprds[[paste0(hicds)]])
    exprds_list
  } # end-foreach iterating over hicds
  names(all_go_ic_list) <- all_hicds
  outFile <- file.path(outFolder, paste0("all_go_ic_list.Rdata"))
  save(all_go_ic_list, file = outFile, version=2)
  #stopifnot(length(unlist(all_go_ic_list, recursive=FALSE)) == length(all_datasets)) => not true because of null ?
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, paste0( "all_go_ic_list.Rdata"))
  cat("... load data\n")
  all_go_ic_list <- get(load(outFile))
}

all_nGO_signif_tads <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_tads"]])))
all_nGO_signif_limma <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_limma"]])))
all_nGO_signif_tadsOnly <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_tadsOnly"]])))
all_nGO_signif_limmaOnly <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_limmaOnly"]])))
all_nGO_signif_intersect <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_intersect"]])))
sum(all_nGO_signif_tads)
sum(all_nGO_signif_limma)
sum(na.omit(all_nGO_signif_tadsOnly))
sum(na.omit(all_nGO_signif_limmaOnly))
sum(na.omit(all_nGO_signif_intersect))
# 
# > sum(all_nGO_signif_tads)
# [1] 3677
# > sum(all_nGO_signif_limma)
# [1] 1870
# > sum(na.omit(all_nGO_signif_tadsOnly))
# [1] 1638
# > sum(na.omit(all_nGO_signif_limmaOnly))
# [1] 1329
# > sum(na.omit(all_nGO_signif_intersect))
# [1] 2800
# > outFolder
# [1] "GO_SPECIFICITY_GENELEVEL_TADLEVEL_INTERSECTDIFF/tadPvalThresh0.01_genePvalThresh0.01"
# > outFolder="GO_SPECIFICITY_GENELEVEL_TADLEVEL_INTERSECTDIFF/tadPvalThresh0.05_genePvalThresh0.05"
# > outFile <- file.path(outFolder, paste0( "all_go_ic_list.Rdata"))
# >   cat("... load data\n")
# ... load data
# >   all_go_ic_list <- get(load(outFile))
# > all_nGO_signif_tads <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_tads"]])))
# > all_nGO_signif_limma <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_limma"]])))
# > all_nGO_signif_tadsOnly <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_tadsOnly"]])))
# > all_nGO_signif_limmaOnly <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_limmaOnly"]])))
# > all_nGO_signif_intersect <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_intersect"]])))
# > sum(all_nGO_signif_tads)
# [1] 2105
# > sum(all_nGO_signif_limma)
# [1] 1294
# > sum(na.omit(all_nGO_signif_tadsOnly))
# [1] 1020
# > sum(na.omit(all_nGO_signif_limmaOnly))
# [1] 508
# > sum(na.omit(all_nGO_signif_intersect))
# [1] 1440
# > 


all_ic_signif_tads <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["signif_tad_go_ic"]])))
all_ic_signif_limma <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["signif_limma_go_ic"]])))
all_ic_signif_tadsOnly <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["signif_tadOnly_go_ic"]])))
all_ic_signif_limmaOnly <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["signif_limmaOnly_go_ic"]])))
all_ic_signif_intersect <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["signif_intersect_go_ic"]])))



######################################################################################## PLOT IC

plot_dt <- rbind(
  data.frame(
    signif_type="tad_signif",
    go_ic = as.numeric(all_ic_signif_tads),
    stringsAsFactors = FALSE
  ),
  data.frame(
    signif_type="limma_signif",
    go_ic = as.numeric(all_ic_signif_limma),
    stringsAsFactors = FALSE
  ),
  data.frame(
    signif_type="tadOnly_signif",
    go_ic = as.numeric(all_ic_signif_tadsOnly),
    stringsAsFactors = FALSE
  ),
  data.frame(
    signif_type="limmaOnly_signif",
    go_ic = as.numeric(all_ic_signif_limmaOnly),
    stringsAsFactors = FALSE
  ),
  data.frame(
    signif_type="intersect",
    go_ic = as.numeric(all_ic_signif_intersect),
    stringsAsFactors = FALSE
  )
)

myLevels <- c("tad_signif", "tadOnly_signif", "limma_signif", "limmaOnly_signif", "intersect")
plot_dt$signif_type <- factor(plot_dt$signif_type, levels=myLevels)
stopifnot(!is.na(plot_dt$signif_type))


# save(plot_dt, file="plot_dt1.Rdata", version=2)

notNa_tads <- sum(!is.na(plot_dt$go_ic[plot_dt$signif_type=="tad_signif"]))
notNa_limma <- sum(!is.na(plot_dt$go_ic[plot_dt$signif_type=="limma_signif"]))
notNa_tadsOnly <- sum(!is.na(plot_dt$go_ic[plot_dt$signif_type=="tadOnly_signif"]))
notNa_limmaOnly <- sum(!is.na(plot_dt$go_ic[plot_dt$signif_type=="limmaOnly_signif"]))
notNa_intersect <- sum(!is.na(plot_dt$go_ic[plot_dt$signif_type=="intersect"]))

subTit <- paste0("notNa_tads=", notNa_tads, "; notNa_limma=", notNa_limma, "; notNa_tadsOnly=",
                 notNa_tadsOnly, "; notNa_limmaOnly=", notNa_limmaOnly, "; notNa_intersect=", notNa_intersect)


p <- ggdensity(plot_dt, 
               title = paste0("enriched GO IC"),
               subtitle=paste0(subTit),
               x = "go_ic", 
               xlab="enriched GO IC",
               color = "signif_type", fill = "signif_type",
               add = "mean", rug = TRUE,
               palette = c(col1,col2,col3,col4,col5))
p <- p+theme(plot.subtitle = element_text(size=10))

outFile <- file.path(outFolder, paste0("all_ds_enrichedGO_IC_", file_suffix, "_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))


p <- ggviolin(plot_dt, 
              title = paste0("enriched GO IC"),
              subtitle=paste0(subTit),
              x="signif_type",
              y = "go_ic", 
              ylab="enriched GO IC",
              xlab="",
              color = "signif_type", 
              # fill = "signif_type",
              add = "mean", rug = TRUE,
              palette = c(col1,col2,col3,col4,col5))
p <- p + scale_y_continuous(breaks = seq(0, max(na.omit(plot_dt$go_ic)), by = 20))
p <- p+theme(plot.subtitle = element_text(size=10))

outFile <- file.path(outFolder, paste0("all_ds_enrichedGO_IC_", file_suffix, "_violin.", plotType))
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
  ),
  data.frame(
    signif_type="tadOnly_signif",
    nbr_GO = as.numeric(all_nGO_signif_tadsOnly),
    stringsAsFactors = FALSE
  ),
  data.frame(
    signif_type="limmaOnly_signif",
    nbr_GO = as.numeric(all_nGO_signif_limmaOnly),
    stringsAsFactors = FALSE
  ),
  data.frame(
    signif_type="intersect",
    nbr_GO = as.numeric(all_nGO_signif_intersect),
    stringsAsFactors = FALSE
  )
)




# save(plot_dt, file="plot_dt2.Rdata", version=2)

myLevels <- c("tad_signif", "tadOnly_signif", "limma_signif", "limmaOnly_signif", "intersect")
plot_dt$signif_type <- factor(plot_dt$signif_type, levels=myLevels)
stopifnot(!is.na(plot_dt$signif_type))

notNa_tads <- sum(!is.na(plot_dt$nbr_GO[plot_dt$signif_type=="tad_signif"]))
notNa_limma <- sum(!is.na(plot_dt$nbr_GO[plot_dt$signif_type=="limma_signif"]))
notNa_tadsOnly <- sum(!is.na(plot_dt$nbr_GO[plot_dt$signif_type=="tadOnly_signif"]))
notNa_limmaOnly <- sum(!is.na(plot_dt$nbr_GO[plot_dt$signif_type=="limmaOnly_signif"]))
notNa_intersect <- sum(!is.na(plot_dt$nbr_GO[plot_dt$signif_type=="intersect"]))

subTit <- paste0("notNa_tads=", notNa_tads, "; notNa_limma=", notNa_limma, "; notNa_tadsOnly=",
                 notNa_tadsOnly, "; notNa_limmaOnly=", notNa_limmaOnly, "; notNa_intersect=", notNa_intersect)

p <- ggdensity(plot_dt, 
               title = paste0("# enriched GO"),
               subtitle=paste0(subTit),
               x = "nbr_GO", 
               xlab="# enriched GO",
               color = "signif_type", fill = "signif_type",
               add = "mean", rug = TRUE,
               palette = c(col1,col2,col3,col4,col5))
p <- p+theme(plot.subtitle = element_text(size=10))

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
              color = "signif_type", 
              # fill = "signif_type",
              add = "mean", rug = TRUE,
              palette = c(col1,col2,col3,col4,col5))
p <- p + scale_y_continuous(breaks = seq(0, max(na.omit(plot_dt$nbr_GO)), by = 20))
p <- p+theme(plot.subtitle = element_text(size=10))

outFile <- file.path(outFolder, paste0("all_ds_nbr_enrichedGO_", file_suffix, "_violin.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))

######################################################################################
stopifnot(names(all_nGO_signif_tads) == names(all_nGO_signif_limma))
stopifnot(names(all_nGO_signif_tads) == names(all_nGO_signif_limmaOnly))
stopifnot(names(all_nGO_signif_tads) == names(all_nGO_signif_tadsOnly))
stopifnot(names(all_nGO_signif_tads) == names(all_nGO_signif_intersect))

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

plot_dt <- as.data.frame(cbind(all_nGO_signif_tads, all_nGO_signif_limma, all_nGO_signif_tadsOnly, all_nGO_signif_limmaOnly, all_nGO_signif_intersect))
plot_dt$exprds <- basename(rownames(plot_dt))

plot_dt$hicds <- gsub("\\.", "/", rownames(plot_dt))

plot_dt$hicds <- basename(dirname(plot_dt$hicds))

plot_dt$ds_label <- paste0(plot_dt$hicds, "\n", plot_dt$exprds)

stopifnot(plot_dt$exprds %in% names(all_cmps))

plot_dt$cmpType <- all_cmps[paste0(plot_dt$exprds)]
# plot_dt$cmpType <- factor(plot_dt$cmpType)
plot_dt$cmpTypeCol <- ifelse(plot_dt$cmpType == "wt_vs_mut", "red", 
                             ifelse(plot_dt$cmpType == "norm_vs_tumor", "blue",
                                    ifelse(plot_dt$cmpType == "subtypes", "green", NA )))
stopifnot(!is.na(plot_dt$cmpTypeCol))

plot_dt$all_nGO_signif_tads_log10 <- log10(plot_dt$all_nGO_signif_tads)
plot_dt$all_nGO_signif_tadsOnly_log10 <- log10(plot_dt$all_nGO_signif_tadsOnly)


plot_dt$all_nGO_signif_limma_log10 <- log10(plot_dt$all_nGO_signif_limma)
plot_dt$all_nGO_signif_limmaOnly_log10 <- log10(plot_dt$all_nGO_signif_limmaOnly)

plot_dt$all_nGO_signif_intersect_log10 <- log10(plot_dt$all_nGO_signif_intersect)


all_x <- c("all_nGO_signif_tads", "all_nGO_signif_tadsOnly", "all_nGO_signif_tads_log10", "all_nGO_signif_tadsOnly_log10")
all_y <- c("all_nGO_signif_limma", "all_nGO_signif_limmaOnly", "all_nGO_signif_limma_log10", "all_nGO_signif_limmaOnly_log10")
myx=all_x[1]
myy=all_y[1]
for(myx in all_x) {
  for(myy in all_y) {
    # p <- ggscatter(plot_dt, x =myx, y = myy)
    p <- ggscatter(plot_dt, x =myx, y = myy, color = "cmpType", label="ds_label", palette=c("wt_vs_mut" = "red", "norm_vs_tumor"="blue", "subtypes"="green"))
    p <- p + rremove("legend")
    # p
    outFile <- file.path(outFolder, paste0("all_ds_", myy, "_vs_", myx, "_scatterplot.", plotType))
    ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile,"\n"))
  }
}

all_x <- c(all_x, all_y)
all_y <- c("all_nGO_signif_intersect", "all_nGO_signif_intersect_log10")

for(myx in all_x) {
  for(myy in all_y) {
    # p <- ggscatter(plot_dt, x =myx, y = myy)
    p <- ggscatter(plot_dt, x =myx, y = myy, color = "cmpType", label="ds_label", palette=c("wt_vs_mut" = "red", "norm_vs_tumor"="blue", "subtypes"="green"))
    p <- p + rremove("legend")
    outFile <- file.path(outFolder, paste0("all_ds_", myy, "_vs_", myx, "_scatterplot.", plotType))
    ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile,"\n"))
  }
}




######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



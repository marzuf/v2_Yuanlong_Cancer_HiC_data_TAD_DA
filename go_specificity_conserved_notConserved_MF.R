script_name <- "go_specificity_conserved_notConserved_MF.R"
options(scipen=100)

SSHFS=F

cat("> START ", script_name, "\n")

# Rscript go_specificity_conserved_notConserved_MF.R 
# Rscript go_specificity_conserved_notConserved_MF.R norm_vs_tumor
# Rscript go_specificity_conserved_notConserved_MF.R subtypes
# Rscript go_specificity_conserved_notConserved_MF.R wt_vs_mut

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

cmpType <- "norm_vs_tumor"
#tads_signifThresh <- args[1]
#genes_signifTresh <- args[2]
if(length(args) == 1) {
  cmpType <- args[1]
} else {
  cmpType <- ""
}


signif_column <- "adjPvalComb"
signifThresh <- 0.01
minOverlapBpRatio <- 0.8
minIntersectGenes <- 3


file_suffix <- paste0(signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes)

outFolder <- file.path("GO_SPECIFICITY_CONSERVED_NOTCONSERVED_MF", cmpType, file_suffix)
dir.create(outFolder, recursive = TRUE)

logFile <- file.path(outFolder, "go_specificity_conserved_notConserved_logFile.txt")
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


inFile <- file.path("GO_SIGNIF_CONSERVED_NOTCONSERVED_MF", cmpType, file_suffix, "all_go_enrich_list.Rdata")
stopifnot(file.exists(inFile))
all_go_enrich_list <- get(load(inFile))

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

if(buildTable) {
  
  hicds = all_hicds[1]
  all_go_ic_list <- foreach(hicds = all_hicds) %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_list <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      cat(paste0("... start: ", hicds, " - ", exprds, "\n"))
      
      go_signif_conserved_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["conserved_signif_tads_genes_resultDT"]]
      
      if(!is.null(go_signif_conserved_dt)) {
        txt <- paste0(hicds, " - ", exprds, " - signif conserved.: # annot. GO:\t", nrow(go_signif_conserved_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_conserved_dt <- go_signif_conserved_dt[go_signif_conserved_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - signif conserved.: # signif. annot. GO:\t", nrow(go_signif_conserved_dt), "\n")
        printAndLog(txt, logFile)
        
        nEnrichedGO_signif_conserved <- nrow(go_signif_conserved_dt)
        
        signif_conserved_go_terms <- rownames(go_signif_conserved_dt)
  
        signif_conserved_go_ids <- my_GO2TERM$go_id[my_GO2TERM$term %in% signif_conserved_go_terms]
        txt <- paste0(hicds, " - ", exprds, " - ... signif conserved.: found GO ids matching GO terms:\t", length(signif_conserved_go_ids), "/", length(signif_conserved_go_terms), "\n")
        printAndLog(txt, logFile)
        
        signif_conserved_go_ic <- GO_IC[names(GO_IC) %in% signif_conserved_go_ids]
        txt <- paste0(hicds, " - ", exprds, " - ... signif conserved.: found GO information content matching GO ids:\t", length(signif_conserved_go_ic), "/", length(signif_conserved_go_ids), "\n")
        printAndLog(txt, logFile)
        txt <- paste0(hicds, " - ", exprds, " - ... signif conserved.: total retrieved data:\t", length(signif_conserved_go_terms), " -> ", length(signif_conserved_go_ic), "\n")
        printAndLog(txt, logFile)
        
        
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... signif conserved.: NULL \n")
        printAndLog(txt, logFile)
        
        nEnrichedGO_signif_conserved <- 0
        signif_conserved_go_ic <- NULL
      }
      go_signif_not_conserved_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["not_conserved_signif_tads_genes_resultDT"]]
      
      if(!is.null(go_signif_not_conserved_dt)) {
        txt <- paste0(hicds, " - ", exprds, " - signif not conserved: # annot. GO:\t", nrow(go_signif_not_conserved_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_not_conserved_dt <- go_signif_not_conserved_dt[go_signif_not_conserved_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - signif not conserved: # signif. annot. GO:\t", nrow(go_signif_not_conserved_dt), "\n")
        printAndLog(txt, logFile)
        
        nEnrichedGO_signif_not_conserved <- nrow(go_signif_not_conserved_dt)
        
        signif_not_conserved_go_terms <- rownames(go_signif_not_conserved_dt)
        
        signif_not_conserved_go_ids <- my_GO2TERM$go_id[my_GO2TERM$term %in% signif_not_conserved_go_terms]
        txt <- paste0(hicds, " - ", exprds, " - ... signif not conserved: found GO ids matching GO terms:\t", length(signif_not_conserved_go_ids), "/", length(signif_not_conserved_go_terms), "\n")
        printAndLog(txt, logFile)
        
        signif_not_conserved_go_ic <- GO_IC[names(GO_IC) %in% signif_not_conserved_go_ids]
        txt <- paste0(hicds, " - ", exprds, " - ... signif not conserved: found GO information content matching GO ids:\t", length(signif_not_conserved_go_ic), "/", length(signif_not_conserved_go_ids), "\n")
        printAndLog(txt, logFile)
        txt <- paste0(hicds, " - ", exprds, " - ... signif not conserved: total retrieved data:\t", length(signif_not_conserved_go_terms), " -> ", length(signif_not_conserved_go_ic), "\n")
        printAndLog(txt, logFile)
        
        
        
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... signif not conserved.: NULL \n")
        printAndLog(txt, logFile)
        nEnrichedGO_signif_not_conserved <- 0
        signif_not_conserved_go_ic <- NULL
      }
      list(
        nEnrichedGO_signif_conserved = nEnrichedGO_signif_conserved,
        nEnrichedGO_signif_not_conserved = nEnrichedGO_signif_not_conserved,
        signif_conserved_go_ic = signif_conserved_go_ic,
        signif_not_conserved_go_ic = signif_not_conserved_go_ic
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

all_nGO_signif_not_conserved <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_not_conserved"]])))
all_nGO_signif_conserved <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["nEnrichedGO_signif_conserved"]])))


all_ic_signif_not_conserved <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["signif_not_conserved_go_ic"]])))
all_ic_signif_conserved <- unlist(lapply(all_go_ic_list, function(sublist) lapply(sublist, function(x) x[["signif_conserved_go_ic"]])))



######################################################################################## PLOT IC

plot_dt <- rbind(
  data.frame(
    signif_type="signif_not_conserved",
    go_ic = as.numeric(all_ic_signif_not_conserved),
    stringsAsFactors = FALSE
  ),
  data.frame(
    signif_type="signif_conserved",
  go_ic = as.numeric(all_ic_signif_conserved),
    stringsAsFactors = FALSE
  )
)

# save(plot_dt, file="plot_dt1.Rdata", version=2)

notNa_signif_conserved <- sum(!is.na(plot_dt$go_ic[plot_dt$signif_type=="signif_conserved"]))
notNa_signif_not_conserved <- sum(!is.na(plot_dt$go_ic[plot_dt$signif_type=="signif_not_conserved"]))

subTit <- paste0("notNa_signif_conserved=", notNa_signif_conserved, "; notNa_signif_not_conserved=", notNa_signif_not_conserved)

p <- ggdensity(plot_dt, 
               title = paste0("enriched GO IC"),
               subtitle=paste0(subTit),
               x = "go_ic", 
               xlab="enriched GO IC",
               color = "signif_type", fill = "signif_type",
               add = "mean", rug = TRUE,
               palette = c(limmaCol, tadCol))

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
               palette = c(limmaCol, tadCol))
p <- p + scale_y_continuous(breaks = seq(0, max(na.omit(plot_dt$go_ic)), by = 20))

outFile <- file.path(outFolder, paste0("all_ds_enrichedGO_IC_", file_suffix, "_violin.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))



######################################################################################## PLOT IC


plot_dt <- rbind(
  data.frame(
    signif_type="signif_not_conserved",
    nbr_GO = as.numeric(all_nGO_signif_not_conserved),
    stringsAsFactors = FALSE
  ),
  data.frame(
    signif_type="signif_conserved",
    nbr_GO = as.numeric(all_nGO_signif_conserved),
    stringsAsFactors = FALSE
  )
)

# save(plot_dt, file="plot_dt2.Rdata", version=2)

notNa_signif_conserved <- sum(!is.na(plot_dt$nbr_GO[plot_dt$signif_type=="signif_conserved"]))
notNa_signif_not_conserved <- sum(!is.na(plot_dt$nbr_GO[plot_dt$signif_type=="signif_not_conserved"]))

subTit <- paste0("notNa_signif_conserved=", notNa_signif_conserved, "; notNa_signif_not_conserved=", notNa_signif_not_conserved)

p <- ggdensity(plot_dt, 
               title = paste0("# enriched GO"),
               subtitle=paste0(subTit),
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
               subtitle=paste0(subTit),
              x="signif_type",
              y = "nbr_GO", 
              ylab="# enriched GO",
              xlab="",
              color = "signif_type", 
              # fill = "signif_type",
              add = "mean", rug = TRUE,
              palette = c(limmaCol, tadCol))
p <- p + scale_y_continuous(breaks = seq(0, max(na.omit(plot_dt$nbr_GO)), by = 20))


outFile <- file.path(outFolder, paste0("all_ds_nbr_enrichedGO_", file_suffix, "_violin.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))


######################################################################################
stopifnot(names(all_nGO_signif_not_conserved) == names(all_nGO_signif_conserved))

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

plot_dt <- as.data.frame(cbind(all_nGO_signif_conserved, all_nGO_signif_not_conserved))
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

plot_dt$all_nGO_signif_conserved_log10 <- log10(plot_dt$all_nGO_signif_conserved)
plot_dt$all_nGO_signif_not_conserved_log10 <- log10(plot_dt$all_nGO_signif_not_conserved)




all_x <- c("all_nGO_signif_conserved", "all_nGO_signif_conserved_log10")
all_y <- c("all_nGO_signif_not_conserved", "all_nGO_signif_not_conserved_log10")
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



######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




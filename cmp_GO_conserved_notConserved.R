startTime <- Sys.time()
cat(paste0("> Rscript cmp_GO_conserved_notConserved.R\n"))

# Rscript cmp_GO_conserved_notConserved.R 
# Rscript cmp_GO_conserved_notConserved.R norm_vs_tumor
# Rscript cmp_GO_conserved_notConserved.R subtypes
# Rscript cmp_GO_conserved_notConserved.R wt_vs_mut


options(scipen=100)

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


pipFolder <- file.path(".")

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

hicds="Panc1_rep12_40kb"
exprds="TCGApaad_wt_mutKRAS"

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight


padjVarGO <- "p.adjust" # p.adjust or qvalue ???


args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  cmpType <- args[1]
} else {
  cmpType <- ""
}

setDir <- "/media/electron"
setDir <- ""


signif_column <- "adjPvalComb"
signifThresh <- 0.01
minOverlapBpRatio <- 0.8
minIntersectGenes <- 3
file_suffix <- paste0(signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes)

outFolder <- file.path("CMP_GO_CONSERVED_NOTCONSERVED", cmpType, file_suffix)
dir.create(outFolder, recursive = TRUE)



inFile <- file.path("GO_SIGNIF_CONSERVED_NOTCONSERVED", cmpType, file_suffix, "all_go_enrich_list.Rdata")
stopifnot(file.exists(inFile))
all_go_enrich_list <- get(load(inFile))

signif_go_pvalThresh <- -log10(0.05)

padjVarGO_plotThresh <- 0.05

all_dt=c(
"conserved_signif_tads_genes_resultDT",
"not_conserved_signif_tads_genes_resultDT"
)

dt = all_dt[1]

topCommonBars <- 10

curr_dataset = names(all_go_enrich_list)[1]


###################################################################################### barplot signif count enriched GO across datasets conserved and not conserved


all_types_signifGO <- foreach(dt = all_dt) %do% {
  all_go_categories <- foreach(curr_dataset = names(all_go_enrich_list)) %dopar% {
      curr_dt <- all_go_enrich_list[[paste0(curr_dataset)]][[paste0(dt)]]
    rownames(curr_dt)[curr_dt[,paste0(padjVarGO)] <= padjVarGO_plotThresh]
  } 
  names(all_go_categories) <- names(all_go_enrich_list)
  go_categories <- unlist(all_go_categories)
  go_categories_count <- setNames(as.numeric(table(go_categories)), names(table(go_categories)))
  go_categories_count <- sort(go_categories_count, decreasing = TRUE)
  outFile <- file.path(outFolder,paste0("all_ds_", dt, "_intersect_",padjVarGO, "_barplot", ".", plotType))
  do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
  par(oma=c(10,1,1,1))
  barplot(go_categories_count[1:topCommonBars], las=2, 
          ylab="# of datasets",
          names.arg=gsub("GO_", "", names(go_categories_count[1:topCommonBars])),
          main=paste0(gsub("_resultDT", "", dt), " - intersect across DS"),
          cex.names=0.6
          )
  
  mtext(side=3, text=paste0(padjVarGO, "<=",padjVarGO_plotThresh, " (top ", topCommonBars, ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder,paste0("all_ds_", dt, "_intersect_",padjVarGO, "_textTable", ".", "txt"))
  write.table(data.frame(go=names(go_categories_count),count=go_categories_count,stringsAsFactors = FALSE), 
              file = outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=F, append=F)
  cat(paste0("... written: ", outFile, "\n"))
  all_go_categories 
}
names(all_types_signifGO) <- all_dt

outFile <- file.path(outFolder,"all_types_signifGO.Rdata")
save(all_types_signifGO, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

###################################################################################### the same by cmpType

if(cmpType == "") {
  all_cmp_names <- unique(all_cmps)
  cmp=all_cmp_names[1]
  for(cmp in all_cmp_names){
    
    toKeepDS <- names(all_go_enrich_list)[basename(names(all_go_enrich_list)) %in% names(all_cmps)[all_cmps==cmp]]
    
    for(dt in all_dt) {
      all_go_categories <- foreach(curr_dataset = toKeepDS) %dopar% {
        curr_dt <- all_go_enrich_list[[paste0(curr_dataset)]][[paste0(dt)]]
        rownames(curr_dt)[curr_dt[,paste0(padjVarGO)] <= padjVarGO_plotThresh]
      } 
      go_categories <- unlist(all_go_categories)
      go_categories_count <- setNames(as.numeric(table(go_categories)), names(table(go_categories)))
      go_categories_count <- sort(go_categories_count, decreasing = TRUE)
      outFile <- file.path(outFolder,paste0("all_ds_",cmp, "_", dt, "_intersect_",padjVarGO, "_barplot", ".", plotType))
      do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
      par(oma=c(10,1,1,1))
      barplot(go_categories_count[1:topCommonBars], las=2, 
              names.arg=gsub("GO_", "", names(go_categories_count[1:topCommonBars])),
              ylab="# of datasets",
              main=paste0(gsub("_resultDT", "", dt), " - ", cmp),
              cex.names=0.6
      )
      mtext(side=3, text=paste0(padjVarGO, "<=",padjVarGO_plotThresh, " (top ", topCommonBars, ")"))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      outFile <- file.path(outFolder,paste0("all_ds_",cmp, "_", dt, "_intersect_",padjVarGO, "_textTable", ".", "txt"))
      write.table(data.frame(go=names(go_categories_count),count=go_categories_count,stringsAsFactors = FALSE), 
                  file = outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=F, append=F)
      cat(paste0("... written: ", outFile, "\n"))
      
    }
  }
    
}

###################################################################################### count common and intersect GO

all_ds <- names(all_types_signifGO[[1]])
stopifnot(setequal(names(all_types_signifGO[[1]]), names(all_types_signifGO[[2]])))


all_ds_GOtypes_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  
  
  conserved_signif_go <- all_types_signifGO[["conserved_signif_tads_genes_resultDT"]][[paste0(ds)]]
  not_conserved_signif_go <- all_types_signifGO[["not_conserved_signif_tads_genes_resultDT"]][[paste0(ds)]]
  
  conservedOnlyGO <- setdiff(conserved_signif_go, not_conserved_signif_go)
  notConservedOnlyGO <- setdiff(not_conserved_signif_go, conserved_signif_go)
  intersectGO <- intersect(conserved_signif_go, not_conserved_signif_go)
  
  if(length(intersectGO) == 0) intersectGO <- NA
  if(length(notConservedOnlyGO) == 0) notConservedOnlyGO <- NA
  if(length(conservedOnlyGO) == 0) conservedOnlyGO <- NA
  
  rbind(
    data.frame(
    dataset = ds,
    GO_type = "conservedOnly",
    GO_term = conservedOnlyGO,
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset = ds,
    GO_type = "notConservedOnly",
    GO_term = notConservedOnlyGO,
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset = ds,
    GO_type = "intersectConservedNotConserved",
    GO_term = intersectGO,
    stringsAsFactors = FALSE
  )
  )
  
}
outFile <- file.path(outFolder,"all_ds_GOtypes_dt.Rdata")
save(all_ds_GOtypes_dt, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))



all_ds_GOtypes_dt_countDS <- aggregate(dataset ~ GO_type + GO_term, FUN=function(x) length(unique(x)), data=all_ds_GOtypes_dt)
gotype = unique(all_ds_GOtypes_dt_countDS$GO_type)[1]
for(gotype in unique(all_ds_GOtypes_dt_countDS$GO_type)) {
  
  sub_dt <- all_ds_GOtypes_dt_countDS[all_ds_GOtypes_dt_countDS$GO_type == gotype,]
  stopifnot(!duplicated(sub_dt$GO_term))
  go_categories_count <- setNames(sub_dt$dataset, sub_dt$GO_term)
  
  go_categories_count <- sort(go_categories_count, decreasing = TRUE)
  outFile <- file.path(outFolder,paste0("all_ds_count","_", gotype, "_", padjVarGO, "_barplot", ".", plotType))
  do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
  par(oma=c(10,1,1,1))
  barplot(go_categories_count[1:topCommonBars], las=2, 
          names.arg=gsub("GO_", "", names(go_categories_count[1:topCommonBars])),
          ylab="# of datasets",
          main=paste0(gotype),
          cex.names=0.6
  )
  mtext(side=3, text=paste0(padjVarGO, "<=",padjVarGO_plotThresh, " (top ", topCommonBars, ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

# PLOT THE GO IC
require(clusterProfiler)
my_GO2TERM <- clusterProfiler:::get_GO2TERM_table()
my_GO2TERM$term <- paste0("GO_", toupper(gsub(" ", "_", my_GO2TERM$Term)))
require(ontologySimilarity)
data(GO_IC)

require(ggpubr)

col1 <- get_palette("Dark2", 5)[1]
col2 <- get_palette("Dark2", 5)[2]
col3 <- get_palette("Dark2", 5)[3]
col4 <- get_palette("Dark2", 5)[4]
col5 <- get_palette("Dark2", 5)[5]

GO_term_id <- setNames(my_GO2TERM$term, my_GO2TERM$go_id)

my_GO_IC <- GO_IC[names(GO_IC) %in% my_GO2TERM$go_id]
stopifnot(names(my_GO_IC) %in% names(GO_term_id))
names(my_GO_IC) <- GO_term_id[paste0(names(my_GO_IC))]

all_ds_GOtypes_dt_ic <- all_ds_GOtypes_dt[all_ds_GOtypes_dt$GO_term %in% names(my_GO_IC),]
all_ds_GOtypes_dt_ic$IC <- my_GO_IC[paste0(all_ds_GOtypes_dt_ic$GO_term)]
stopifnot(!is.na(all_ds_GOtypes_dt_ic$IC))

all_ds_GOtypes_dt_ic_count <- aggregate(GO_term ~ dataset+GO_type, data=all_ds_GOtypes_dt_ic, FUN=length)

notNa_conservedOnly <- sum(all_ds_GOtypes_dt_ic$GO_type=="conservedOnly")
notNa_notConservedOnly <- sum(all_ds_GOtypes_dt_ic$GO_type=="notConservedOnly")
notNa_intersectConservedNotConserved <- sum(all_ds_GOtypes_dt_ic$GO_type=="intersectConservedNotConserved")

subTit <- paste0("notNa_conservedOnly=", notNa_conservedOnly, "; notNa_notConservedOnly=", notNa_notConservedOnly, "; notNa_intersect=", notNa_intersectConservedNotConserved)

p <- ggdensity(all_ds_GOtypes_dt_ic, 
               title = paste0("enriched GO IC"),
               subtitle=paste0(subTit),
               x = "IC", 
               xlab="enriched GO IC",
               color = "GO_type", fill = "GO_type",
               # add = "mean", rug = TRUE,
               palette = c(col1, col2, col3))

outFile <- file.path(outFolder, paste0("all_ds_enrichedGO_IC_byGOtype_", file_suffix, "_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))

notNa_conservedOnly <- sum(all_ds_GOtypes_dt_ic_count$GO_type=="conservedOnly")
notNa_notConservedOnly <- sum(all_ds_GOtypes_dt_ic_count$GO_type=="notConservedOnly")
notNa_intersectConservedNotConserved <- sum(all_ds_GOtypes_dt_ic_count$GO_type=="intersectConservedNotConserved")

subTit <- paste0("notNa_conservedOnly=", notNa_conservedOnly, "; notNa_notConservedOnly=", notNa_notConservedOnly, "; notNa_intersect=", notNa_intersectConservedNotConserved)


p <- ggdensity(all_ds_GOtypes_dt_ic_count, 
               title = paste0("# enriched GO IC"),
               subtitle=paste0(subTit),
               x = "GO_term", 
               xlab="# enriched GO IC",
               color = "GO_type", fill = "GO_type",
               # add = "mean", rug = TRUE,
               palette = c(col1, col2, col3))

outFile <- file.path(outFolder, paste0("all_ds_nbrEnrichedGO_IC_byGOtype_", file_suffix, "_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))







############################################################################################################################################################################ same by cmpType

if(cmpType == "") {
  
  all_cmp_names <- unique(all_cmps)
  cmp=all_cmp_names[2]
  for(cmp in all_cmp_names){
    
    toKeepDS <- names(all_go_enrich_list)[basename(names(all_go_enrich_list)) %in% names(all_cmps)[all_cmps==cmp]]
    
    cmp_all_ds_GOtypes_dt <- all_ds_GOtypes_dt[all_ds_GOtypes_dt$dataset %in% toKeepDS,]
    
    all_ds_GOtypes_dt_countDS <- aggregate(dataset ~ GO_type + GO_term, FUN=function(x) length(unique(x)), data=cmp_all_ds_GOtypes_dt)
    
    
    
    for(gotype in unique(all_ds_GOtypes_dt_countDS$GO_type)) {
      
      sub_dt <- all_ds_GOtypes_dt_countDS[all_ds_GOtypes_dt_countDS$GO_type == gotype,]
      stopifnot(!duplicated(sub_dt$GO_term))
      go_categories_count <- setNames(sub_dt$dataset, sub_dt$GO_term)
      
      go_categories_count <- sort(go_categories_count, decreasing = TRUE)
      outFile <- file.path(outFolder,paste0("all_ds_count","_",cmp, "_",  gotype, "_", padjVarGO, "_barplot", ".", plotType))
      do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
      par(oma=c(10,1,1,1))
      barplot(go_categories_count[1:topCommonBars], las=2, 
              names.arg=gsub("GO_", "", names(go_categories_count[1:topCommonBars])),
              ylab="# of datasets",
              main=paste0(gotype, " - ", cmp),
              cex.names=0.6
      )
      mtext(side=3, text=paste0(padjVarGO, "<=",padjVarGO_plotThresh, " (top ", topCommonBars, ")"))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      
      
      
    }
    
    
    cmp_all_ds_GOtypes_dt_ic <- cmp_all_ds_GOtypes_dt[cmp_all_ds_GOtypes_dt$GO_term %in% names(my_GO_IC),]
    cmp_all_ds_GOtypes_dt_ic$IC <- my_GO_IC[paste0(cmp_all_ds_GOtypes_dt_ic$GO_term)]
    stopifnot(!is.na(cmp_all_ds_GOtypes_dt_ic$IC))
    
    cmp_all_ds_GOtypes_dt_ic_count <- aggregate(GO_term ~ dataset+GO_type, data=cmp_all_ds_GOtypes_dt_ic, FUN=length)
    
    notNa_conservedOnly <- sum(cmp_all_ds_GOtypes_dt_ic$GO_type=="conservedOnly")
    notNa_notConservedOnly <- sum(cmp_all_ds_GOtypes_dt_ic$GO_type=="notConservedOnly")
    notNa_intersectConservedNotConserved <- sum(cmp_all_ds_GOtypes_dt_ic$GO_type=="intersectConservedNotConserved")
    
    subTit <- paste0("notNa_conservedOnly=", notNa_conservedOnly, "; notNa_notConservedOnly=", notNa_notConservedOnly, "; notNa_intersect=", notNa_intersectConservedNotConserved)
    
    p <- ggdensity(cmp_all_ds_GOtypes_dt_ic, 
                   title = paste0("enriched GO IC - ", cmp),
                   subtitle=paste0(subTit),
                   x = "IC", 
                   xlab="enriched GO IC",
                   color = "GO_type", fill = "GO_type",
                   # add = "mean", rug = TRUE,
                   palette = c(col1, col2, col3))
    
    outFile <- file.path(outFolder, paste0(cmp, "_enrichedGO_IC_byGOtype_", file_suffix, "_density.", plotType))
    ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile,"\n"))
    
    notNa_conservedOnly <- sum(cmp_all_ds_GOtypes_dt_ic_count$GO_type=="conservedOnly")
    notNa_notConservedOnly <- sum(cmp_all_ds_GOtypes_dt_ic_count$GO_type=="notConservedOnly")
    notNa_intersectConservedNotConserved <- sum(cmp_all_ds_GOtypes_dt_ic_count$GO_type=="intersectConservedNotConserved")
    
    subTit <- paste0("notNa_conservedOnly=", notNa_conservedOnly, "; notNa_notConservedOnly=", notNa_notConservedOnly, "; notNa_intersect=", notNa_intersectConservedNotConserved)
    
    
    p <- ggdensity(cmp_all_ds_GOtypes_dt_ic_count, 
                   title = paste0("# enriched GO IC - ", cmp),
                   subtitle=paste0(subTit),
                   x = "GO_term", 
                   xlab="# enriched GO IC",
                   color = "GO_type", fill = "GO_type",
                   # add = "mean", rug = TRUE,
                   palette = c(col1, col2, col3))
    
    outFile <- file.path(outFolder, paste0(cmp, "_nbrEnrichedGO_IC_byGOtype_", file_suffix, "_density.", plotType))
    ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile,"\n"))
    
    
    
    
    
  }
  
}


######################################################################################
######################################################################################
######################################################################################

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


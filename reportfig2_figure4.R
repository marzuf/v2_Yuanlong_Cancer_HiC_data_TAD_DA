startTime <- Sys.time()
cat(paste0("> Rscript reportfig2_figure4.R\n"))

# Rscript reportfig2_figure4.R 0.01 0.05

options(scipen=100)

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


outFolder <- file.path("REPORTFIG2_FIGURE4")
dir.create(outFolder, recursive = TRUE)

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
axisCex <- 1.4

padjVarGO <- "p.adjust" # p.adjust or qvalue ???



TAD_pvalThresh <- 0.01
gene_pvalThresh <- 0.05
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
TAD_pvalThresh <- args[1]
gene_pvalThresh <- args[2]
if(length(args) == 3) {
  cmpType <- args[3]  
} else {
  cmpType <- ""
}



setDir <- "/media/electron"
setDir <- ""


inFile <- file.path("GO_SIGNIF_GENELEVEL_TADLEVEL_INTERSECTDIFF", paste0("tadPvalThresh", TAD_pvalThresh, "_genePvalThresh", gene_pvalThresh), "all_go_enrich_list.Rdata")
all_go_enrich_list <- get(load(inFile))


signif_go_pvalThresh <- -log10(0.05)

padjVarGO_plotThresh <- 0.05

all_dt=c(
  "tad_signif_enrich_resultDT",
  "limma_signif_enrich_resultDT"
)

dt = all_dt[1]

topCommonBars <- 10

curr_dataset = names(all_go_enrich_list)[1]

############################################################################################################################################################################ barplot count GO across ds

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
          cex.lab=axisCex,
          cex.axis=axisCex,
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


############################################################################################################################################################################ same by cmpType

if(cmpType == "") {
  
  all_cmp_names <- unique(all_cmps)
  cmp=all_cmp_names[2]
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
              cex.names=0.6,
              cex.axis=axisCex,
              cex.lab = axisCex
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

# "tadsOnly_signif_enrich_resultDT",
# "limmaOnly_signif_enrich_resultDT",
# "tad_signif_enrich_resultDT",
# "limma_signif_enrich_resultDT",
# "intersect_signif_enrich_resultDT"



all_ds_GOtypes_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  
  tad_signif_go <- all_types_signifGO[["tad_signif_enrich_resultDT"]][[paste0(ds)]]
  limma_signif_go <- all_types_signifGO[["limma_signif_enrich_resultDT"]][[paste0(ds)]]
  
  if(length(tad_signif_go) == 0) tad_signif_go <- NA
  if(length(limma_signif_go) == 0) limma_signif_go <- NA
  
  rbind(
    data.frame(
      dataset = ds,
      GO_type = "tad_signif_GO",
      GO_term = tad_signif_go,
      stringsAsFactors = FALSE
    ),
    data.frame(
      dataset = ds,
      GO_type = "limma_signif_GO",
      GO_term = limma_signif_go,
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
          cex.names=0.6,
          cex.axis=axisCex,
          cex.lab=axisCex
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

col1 <- get_palette("Dark2", 6)[1]
col2 <- get_palette("Dark2", 6)[2]


GO_term_id <- setNames(my_GO2TERM$term, my_GO2TERM$go_id)

my_GO_IC <- GO_IC[names(GO_IC) %in% my_GO2TERM$go_id]
stopifnot(names(my_GO_IC) %in% names(GO_term_id))
names(my_GO_IC) <- GO_term_id[paste0(names(my_GO_IC))]

all_ds_GOtypes_dt_ic <- all_ds_GOtypes_dt[all_ds_GOtypes_dt$GO_term %in% names(my_GO_IC),]
all_ds_GOtypes_dt_ic$IC <- my_GO_IC[paste0(all_ds_GOtypes_dt_ic$GO_term)]
stopifnot(!is.na(all_ds_GOtypes_dt_ic$IC))

all_ds_GOtypes_dt_ic_count <- aggregate(GO_term ~ dataset+GO_type, data=all_ds_GOtypes_dt_ic, FUN=length)


notNa_tad <- sum(all_ds_GOtypes_dt_ic$GO_type=="tad_signif_GO")
notNa_limma <- sum(all_ds_GOtypes_dt_ic$GO_type=="limma_signif_GO")

subTit <- paste0("# signif. TAD level = ", notNa_tad, "\n# signif. gene level = ", notNa_limma)


all_ds_GOtypes_dt_ic$GO_type[all_ds_GOtypes_dt_ic$GO_type == "tad_signif_GO"] <- "gene level signif."
all_ds_GOtypes_dt_ic$GO_type[all_ds_GOtypes_dt_ic$GO_type == "limma_signif_GO"] <- "TAD level signif."

p <- ggdensity(all_ds_GOtypes_dt_ic, 
               title = paste0("Information content of enriched GO"),
               subtitle=paste0(subTit),
               x = "IC", 
               # xlab="enriched GO IC",
               xlab="IC of signif. enriched GO",
               color = "GO_type", fill = "GO_type",
               # add = "mean", rug = TRUE,
               palette = c(col1, col2))+ 
  labs(fill="", col="")+
  theme(plot.title = element_text(hjust=0.5, size=16, face="bold"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14)
  )

outFile <- file.path(outFolder, paste0("all_ds_enrichedGO_IC_byGOtype_density.", "png"))
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
              cex.names=0.6,
              cex.axis=axisCex,
              cex.lab=axisCex
      )
      mtext(side=3, text=paste0(padjVarGO, "<=",padjVarGO_plotThresh, " (top ", topCommonBars, ")"))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
    }
    
    cmp_all_ds_GOtypes_dt_ic <- cmp_all_ds_GOtypes_dt[cmp_all_ds_GOtypes_dt$GO_term %in% names(my_GO_IC),]
    cmp_all_ds_GOtypes_dt_ic$IC <- my_GO_IC[paste0(cmp_all_ds_GOtypes_dt_ic$GO_term)]
    stopifnot(!is.na(cmp_all_ds_GOtypes_dt_ic$IC))
    
    cmp_all_ds_GOtypes_dt_ic_count <- aggregate(GO_term ~ dataset+GO_type, data=cmp_all_ds_GOtypes_dt_ic, FUN=length)
    
    
    notNa_tad <- sum(cmp_all_ds_GOtypes_dt_ic$GO_type=="tad_signif_GO")
    notNa_limma <- sum(cmp_all_ds_GOtypes_dt_ic$GO_type=="limma_signif_GO")
    
    subTit <- paste0("# signif. TAD level = ", notNa_tad, "\n# signif. gene level = ", notNa_limma)
    
    cmp_all_ds_GOtypes_dt_ic$GO_type[cmp_all_ds_GOtypes_dt_ic$GO_type == "tad_signif_GO"] <- "gene level signif."
    cmp_all_ds_GOtypes_dt_ic$GO_type[cmp_all_ds_GOtypes_dt_ic$GO_type == "limma_signif_GO"] <- "TAD level signif."
    
    
    p <- ggdensity(cmp_all_ds_GOtypes_dt_ic, 
                   # title = paste0("enriched GO IC - ", cmp),
                   title = paste0("Information content of enriched GO - ", cmp),
                   subtitle=paste0(subTit),
                   x = "IC", 
                   # xlab="enriched GO IC",
                   xlab="IC of signif. enriched GO",
                   color = "GO_type", fill = "GO_type",
                   # add = "mean", rug = TRUE,
                   palette = c(col1, col2)) + 
      labs(fill="", col="")+
      theme(plot.title = element_text(hjust=0.5, size=16, face="bold"),
            axis.text = element_text(size=12),
            axis.title = element_text(size=14)
            )
    outFile <- file.path(outFolder, paste0(cmp, "_enrichedGO_IC_byGOtype_density.", plotType))
    ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile,"\n"))
  }
}


######################################################################################
######################################################################################
######################################################################################

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

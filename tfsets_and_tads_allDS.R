startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)
  
# Rscript tfsets_and_tads_allDS.R crisp 
# Rscript tfsets_and_tads_allDS.R crisp LG1_40kb
# Rscript tfsets_and_tads_allDS.R c3.mir
# Rscript tfsets_and_tads_allDS.R c3.tft
# Rscript tfsets_and_tads_allDS.R c3.all
# Rscript tfsets_and_tads_allDS.R trrust
# Rscript tfsets_and_tads_allDS.R tftg
# Rscript tfsets_and_tads_allDS.R motifmap LG1_40kb
# Rscript tfsets_and_tads_allDS.R motifmap 
# Rscript tfsets_and_tads_allDS.R kegg LG1_40kb

plotType <- "png"
myHeight <- 400
myWidth <- 600
plotCex <- 1.4


dsIn <- "crisp"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1 | length(args) == 2)
dsIn <- args[1]
if(length(args) == 2) {
  all_hicds <- args[2]
} else {
  all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
}

stopifnot(dsIn %in% c("crisp", "c3.mir", "c3.all", "c3.tft", "trrust", "tftg", "motifmap", "kegg"))

nPermut <- 10000

outFolder <- file.path(paste0("TFSETS_AND_TADS_", toupper(dsIn)))

buildData <- FALSE

for(hicds in all_hicds){
  
  cat(paste0("> START - ", hicds, "\n"))
    
  if(buildData) {
    
    setDir <- "/media/electron"
    setDir <- ""
    entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
    gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
    gff_dt$entrezID <- as.character(gff_dt$entrezID)
    stopifnot(!duplicated(gff_dt$entrezID))
    stopifnot(!duplicated(gff_dt$symbol))
    entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
    symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)
    
    
    if(dsIn == "crisp") {
      reg_file <- file.path("gene_set_library_crisp_processed.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
      cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
      reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
    } else if(dsIn == "trrust"){
      reg_file <- file.path("trrust_rawdata.human.tsv")
      reg_dt <- read.delim(reg_file, sep="\t", header=FALSE, stringsAsFactors = FALSE,
                           col.names = c("regSymbol", "targetSymbol", "direction", "ID"))
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
      cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
      reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
    } else if(dsIn == "tftg") {
      reg_file <- file.path("tftg_db_all_processed.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE, 
                           col.names=c("regSymbol", "targetEntrezID"))
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    } else if(dsIn == "motifmap"){
      reg_file <- file.path("MOTIFMAP_ALLGENES/overlapDT_bp.Rdata")
      reg_dt <- get(load(reg_file))
      colnames(reg_dt)[colnames(reg_dt)=="entrezID"] <- "targetEntrezID"
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    } else if(dsIn == "kegg"){
      
      
      reg_file <- file.path("hsa_kegg_entrez.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=FALSE, stringsAsFactors = FALSE,
                           col.names = c("targetEntrezID", "regSymbol"))
      reg_dt$targetEntrezID <- gsub("hsa:", "",reg_dt$targetEntrezID )
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))

      
      
    }else {
      reg_file <- file.path(paste0(dsIn, ".v7.0.entrez_processed.txt"))
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE, 
                           col.names=c("regSymbol", "targetEntrezID"))
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    }
    
    g2t_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    
    g2t_dt <- read.delim(g2t_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    
    g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    stopifnot(!duplicated(g2t_dt$entrezID))
    g2t_vect <- setNames(g2t_dt$region, g2t_dt$entrezID)
    
    reg_dt <- reg_dt[reg_dt$targetEntrezID %in% g2t_dt$entrezID,]
    cat(paste0("with g2t assignment: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    reg_dt$targetRegion <- g2t_vect[paste0(reg_dt$targetEntrezID)]
    stopifnot(!is.na(reg_dt))
    
    reg_nGenes <- unlist(by(reg_dt, reg_dt$regSymbol, nrow))
    stopifnot(sum(reg_nGenes) == nrow(reg_dt))
    
    reg_nTADs <- unlist(by(reg_dt, reg_dt$regSymbol, function(x) length(unique(x$targetRegion))))
    
    all_gene_nbrs <- unique(reg_nGenes)
    
    stopifnot(length(all_gene_nbrs) > 1)
    
    ng=all_gene_nbrs[1]
    permut_dt <- foreach(i_ng = 1:length(all_gene_nbrs), .combine='rbind') %do% {
      ng <-  all_gene_nbrs[i_ng]
      ng_permut <- foreach(npermut = 1:nPermut, .combine='c') %dopar% {
        cat(paste0(i_ng, "/", length(all_gene_nbrs), " - ", npermut, "/", nPermut, "\n"))
        permut_rows <- sample(1:nrow(reg_dt), size=ng)
        length(unique(reg_dt$targetRegion[permut_rows]))
      }
      data.frame(
        reg_nGenes = ng,
        reg_nTADs_permut = ng_permut,
        stringsAsFactors = FALSE
      )
    }
    
    stopifnot(names(reg_nGenes) == names(reg_nTADs))
    
    obs_dt <- data.frame(
      reg_nGenes = as.numeric(reg_nGenes),
      reg_nTADs_obs = as.numeric(reg_nTADs), 
      stringsAsFactors = FALSE
    )
    outFile <- file.path(outFolder, hicds, "permut_dt.Rdata")
    dir.create(dirname(outFile), recursive = TRUE)
    save(permut_dt, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFolder, hicds, "obs_dt.Rdata")
    dir.create(dirname(outFile), recursive = TRUE)
    save(obs_dt, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
  } else {
    inFile <- file.path(outFolder, hicds, "permut_dt.Rdata")
    permut_dt <- get(load(inFile))
    inFile <- file.path(outFolder, hicds, "obs_dt.Rdata")
    obs_dt <- get(load(inFile))
  }
  
  
  permut_dt <- permut_dt[order(permut_dt$reg_nGenes),]
  
  yMin <- min(c(permut_dt$reg_nTADs_permut, obs_dt$reg_nTADs), na.rm=TRUE)
  yMax <- max(c(permut_dt$reg_nTADs_permut, obs_dt$reg_nTADs), na.rm=TRUE)
  
  outFile <- file.path(outFolder, hicds, paste0("obs_permut_boxplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(reg_nTADs_permut ~ reg_nGenes, data = permut_dt, 
          xlab = "# regulated genes ",
          ylab = "# of TADs",
          # ylim = c(yMin, yMax),
          ylim = c(0, yMax),
          main = paste0(hicds), 
          at = unique(permut_dt$reg_nGenes))
  points(
    x = obs_dt$reg_nGenes,
    y = obs_dt$reg_nTADs,
    cex.lab = plotCex,
    cex.axis = plotCex,
    pch = 16,
    col ="red"
  )
  legend("topleft", 
         legend=c("obs.", "permut."),
         col = c("red", "black"),
         pch=16,
         bty="n")
  mtext(side=3, text = paste0(dsIn))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  permut_dt$tad_ratio <- permut_dt$reg_nTADs_permut/permut_dt$reg_nGenes
  obs_dt$tad_ratio <- obs_dt$reg_nTADs_obs/obs_dt$reg_nGenes
  
  outFile <- file.path(outFolder, hicds, paste0("obs_permut_multidens.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  
  plot_multiDens(
    list(obs = obs_dt$tad_ratio,
         permut = permut_dt$tad_ratio),
    plotTit = paste0(hicds)
  )
  mtext(side=3, text = paste0(dsIn))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  yMin <- min(c(permut_dt$reg_nTADs_permut[permut_dt$reg_nGenes <= 100], obs_dt$reg_nTADs[obs_dt$reg_nGenes<=100]), na.rm=TRUE)
  yMax <- max(c(permut_dt$reg_nTADs_permut[permut_dt$reg_nGenes <= 100], obs_dt$reg_nTADs[obs_dt$reg_nGenes<=100]), na.rm=TRUE)
  
  
  outFile <- file.path(outFolder, hicds, paste0("obs_permut_boxplot_to100.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(reg_nTADs_permut ~ reg_nGenes, data = permut_dt[permut_dt$reg_nGenes <= 100,], 
          xlab = "# regulated genes ",
          ylab = "# of TADs",
          ylim = c(0, yMax),
          main = paste0(hicds), 
          at = unique(permut_dt$reg_nGenes[permut_dt$reg_nGenes <= 100]))
  points(
    x = obs_dt$reg_nGenes[obs_dt$reg_nGenes<=100],
    y = obs_dt$reg_nTADs[obs_dt$reg_nGenes<=100],
    cex.lab = plotCex,
    cex.axis = plotCex,
    pch = 16,
    col ="red"
  )
  mtext(paste0(dsIn, "reg_nGenes <= 100"), side=3)
  legend("topleft", 
         legend=c("obs.", "permut."),
         col = c("red", "black"),
         pch=16,
         bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  permut_dt$tad_ratio <- permut_dt$reg_nTADs_permut/permut_dt$reg_nGenes
  obs_dt$tad_ratio <- obs_dt$reg_nTADs_obs/obs_dt$reg_nGenes
  
  outFile <- file.path(outFolder, hicds, paste0("obs_permut_multidens_to100.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  
  plot_multiDens(
    list(obs = obs_dt$tad_ratio[obs_dt$reg_nGenes <= 100],
         permut = permut_dt$tad_ratio[permut_dt$reg_nGenes <= 100]),
    plotTit = paste0(hicds)
  )
  mtext(paste0(dsIn, "reg_nGenes <= 100"), side=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))

}

#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))

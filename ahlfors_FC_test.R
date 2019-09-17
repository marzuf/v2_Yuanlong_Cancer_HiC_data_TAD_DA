startTime <- Sys.time()
cat(paste0("> Rscript ahlfors_FC_test.R\n"))

plotType <- "png"
myHeight  <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight

# Rscript ahlfors_FC_test.R

script_name <- "ahlfors_FC_test.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS=F
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

buildTable <- TRUE

hicds = "Panc1_rep12_40kb"
exprds = "TCGApaad_wt_mutKRAS"

mainFolder <- file.path(".")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"


randomQtCol <- "lightblue"
nRandom <- 100
set.seed(20190917)
genomeCol <- "black"
genomeTadCol <- "red"

minQt <- 0.05
maxQt <- 0.95

outFolder <- file.path("AHLFORS_FC_TEST")
dir.create(outFolder, recursive = TRUE)

pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))

mainFolder <- file.path(".")
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")

all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

####################################################################################################################################### >>> prepare the data
if(buildTable) {
  all_cumsum_data <- foreach(hicds = all_hicds) %dopar% {


    gene2tad_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(gene2tad_file))
    all_g2t_dt <- read.delim(gene2tad_file, stringsAsFactors = FALSE, col.names = c("entrezID", "chromo", "start", "end", "region"), header=FALSE)
    all_g2t_dt$entrezID <- as.character(all_g2t_dt$entrezID)
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_cumsum <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      geneListFile <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneListFile))
      geneList <- get(load(geneListFile))
      
      topTableFile <- file.path(pipFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")
      stopifnot(file.exists(topTableFile))
      topTable_DT <- get(load(topTableFile))
      topTable_DT$genes <- as.character(topTable_DT$genes) 
      stopifnot(names(geneList) %in%  topTable_DT$genes)
      topTable_DT <- topTable_DT[topTable_DT$genes %in% names(geneList),]
      topTable_DT$entrezID <- geneList[paste0(topTable_DT$genes)]
      stopifnot(!is.na(topTable_DT$entrezID))
      stopifnot(!duplicated(topTable_DT$entrezID))
      
      stopifnot(geneList %in% all_g2t_dt$entrezID)
      g2t_dt <- all_g2t_dt[all_g2t_dt$entrezID %in% geneList,]
      g2t_dt <- g2t_dt[order(g2t_dt$chromo, g2t_dt$start, g2t_dt$end),]
      
      gene_logFC <- setNames(topTable_DT$logFC, topTable_DT$entrezID)
      gene_tad <- setNames(g2t_dt$region, g2t_dt$entrezID)
      stopifnot(setequal(names(gene_logFC), names(gene_tad)))
      
      i=2
      all_sumFC <- foreach(i = 2:c(nrow(g2t_dt)-1)) %dopar% {
        chromo_i <- g2t_dt$chromo[i]
        chromo_next <- g2t_dt$chromo[i+1]
        chromo_prev <- g2t_dt$chromo[i-1]
        
        if( chromo_i != chromo_next | chromo_i != chromo_prev) return(NULL)
        entrez_i <- g2t_dt$entrezID[i]
        stopifnot(entrez_i %in% names(gene_logFC)); stopifnot(entrez_i %in% names(gene_tad))
        entrez_prev <- g2t_dt$entrezID[i-1]
        stopifnot(entrez_prev %in% names(gene_logFC)); stopifnot(entrez_prev %in% names(gene_tad))
        entrez_next <- g2t_dt$entrezID[i-1]
        stopifnot(entrez_next %in% names(gene_logFC)); stopifnot(entrez_next %in% names(gene_tad))
        
        fc_i <- gene_logFC[paste0(entrez_i)]
        fc_prev <- gene_logFC[paste0(entrez_prev)]
        fc_next <- gene_logFC[paste0(entrez_next)]
        
        tad_i <- gene_tad[paste0(entrez_i)]
        tad_prev <- gene_tad[paste0(entrez_prev)]
        tad_next <- gene_tad[paste0(entrez_next)]
        
        sum_fc <- fc_i*fc_prev + fc_i*fc_next
        
        if(tad_i == tad_prev & tad_i == tad_next) {
          i_tad_sum_fc <- sum_fc
        } else {
          i_tad_sum_fc <- NULL
        }
        list(
          i_sum_fc = sum_fc,
          i_tad_sum_fc = i_tad_sum_fc)
      }
      
      all_sumFC_values <- unlist(lapply(all_sumFC, function(x)x[["i_sum_fc"]]))
      cumsum_genome <- cumsum(all_sumFC_values)
      
      all_tad_sumFC_values <- unlist(lapply(all_sumFC, function(x)x[["i_tad_sum_fc"]]))
      cumsum_genome_tads <- cumsum(all_tad_sumFC_values)
      
      nBoundaries <-  length(cumsum_genome) - length(cumsum_genome_tads)
      
      all_random_cumsum <- foreach(irand = 1:nRandom, .combine='rbind') %dopar% {
        random_null_idx <- sample(1:length(cumsum_genome), size=nBoundaries, replace = FALSE)
        random_values <- all_sumFC_values[-random_null_idx]
        stopifnot(length(random_values) == length(all_tad_sumFC_values))
        cumsum_rand <- cumsum(random_values)
        #lines(x=1:length(cumsum_rand), y=cumsum_rand, type="l", col="grey")
        cumsum_rand
      }
      # for each column: the 0.05 and 0.95
      random_qts <- apply(all_random_cumsum, 2, quantile, probs=c(minQt,maxQt))
      
      
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", "cumsumcurves.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      plot(cumsum_genome, type="l", col=genomeCol,
           xlab="genome index",
           ylab="cumsum FC product window",
           main=paste0(hicds, " - ", exprds)
           )
      lines(x=1:length(cumsum_genome_tads), y=cumsum_genome_tads, type="l", col=genomeTadCol)
      polygon( c(1:ncol(random_qts), rev(1:ncol(random_qts))), c(random_qts[1,], rev(random_qts[2,])), col=randomQtCol, 
               density = c(20, 60), angle = c(-45, -45))
      legend("topleft", lty=1, legend=c("no TADs", "with TADs", "with random TADs"), col = c(genomeCol, genomeTadCol, randomQtCol), bty="n", cex=1)
      foo <- dev.off()
      cat(paste0("... written: ", outFile,"\n"))
      
      list(
        cumsum_genome=cumsum_genome,
        cumsum_genome_tads=cumsum_genome_tads,
        random_qts005 = random_qts[1,], 
        random_qts095 = random_qts[2,]
      )
      
      
    } # end-foreach exprds
    exprds_cumsum
  } # end-foreach hicds
  
    
  
  outFile <- file.path(outFolder, paste0("all_cumsum_data.Rdata"))
  save(all_cumsum_data, file=outFile, version=2)  
  cat(paste0("... written: ", outFile, "\n"))

} else {
  # outFile= "AHLFORS_FC_TEST/all_cumsum_data.Rdata"
  outFile <- file.path(outFolder, paste0("all_cumsum_data.Rdata"))
  cat(paste0("... load data\n"))
  all_cumsum_data <- get(load(outFile))
  
}
all_cumsum_data_ul <- unlist(all_cumsum_data, recursive=FALSE)
all_cumsum_genome <- lapply(all_cumsum_data_ul,function(x)x[["cumsum_genome"]])
all_cumsum_genome_tads <- lapply(all_cumsum_data_ul,function(x)x[["cumsum_genome_tads"]])
all_cumsum_random_qts005 <- lapply(all_cumsum_data_ul,function(x)x[["random_qts005"]])
all_cumsum_random_qts095 <- lapply(all_cumsum_data_ul,function(x)x[["random_qts095"]])

minLength <- min(c( lengths(all_cumsum_genome),
                    lengths(all_cumsum_genome_tads),
                    lengths(all_cumsum_random_qts005),
                    lengths(all_cumsum_random_qts095)))

all_cumsum_genome_sub <- do.call(rbind, lapply(all_cumsum_genome, function(x) x[1:minLength]))
all_cumsum_genome_tads_sub <- do.call(rbind, lapply(all_cumsum_genome_tads, function(x) x[1:minLength]))
all_cumsum_random_qts005_sub <- do.call(rbind, lapply(all_cumsum_random_qts005, function(x) x[1:minLength]))
all_cumsum_random_qts095_sub <- do.call(rbind, lapply(all_cumsum_random_qts095, function(x) x[1:minLength]))

maxCumsum <- max(c( unlist(all_cumsum_genome_sub),
                    unlist(all_cumsum_genome_tads_sub),
                    unlist(all_cumsum_random_qts005_sub),
                    unlist(all_cumsum_random_qts095_sub)))

minCumsum <- min(c( unlist(all_cumsum_genome_sub),
                    unlist(all_cumsum_genome_tads_sub),
                    unlist(all_cumsum_random_qts005_sub),
                    unlist(all_cumsum_random_qts095_sub)))



cs_genome_qts <-  apply(all_cumsum_genome_sub, 2, quantile, probs=c(minQt,maxQt))
cs_genome_tads_qts <-  apply(all_cumsum_genome_tads_sub, 2, quantile, probs=c(minQt,maxQt))

cs_genome_mean <-  colMeans(all_cumsum_genome_sub)
cs_genome_tads_mean <-  colMeans(all_cumsum_genome_tads_sub)
# cs_random005_qts <-  apply(all_cumsum_random_qts005_sub, 2, quantile, probs=c(minQt,maxQt))
# cs_random095_qts <-  apply(all_cumsum_random_qts095_sub, 2, quantile, probs=c(minQt,maxQt))
cs_random005_mean <-  colMeans(all_cumsum_random_qts005_sub)
cs_random095_mean <-  colMeans(all_cumsum_random_qts095_sub)
stopifnot(dim(cs_genome_qts) == dim(cs_genome_tads_qts))
# stopifnot(dim(cs_genome_qts) == dim(cs_random005_qts))
# stopifnot(dim(cs_genome_qts) == dim(cs_random095_qts))
stopifnot(ncol(cs_genome_qts) == length(cs_random005_mean))
stopifnot(ncol(cs_genome_qts) == length(cs_random095_mean))

nDS <- nrow(all_cumsum_genome_sub)

# outFile <- file.path(outFolder, paste0("all_ds_cumsumcurves_qts.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# plot(NULL,
# #      type="l", col=genomeCol,
#      xlim=c(0,minLength),
#      ylim=c(minCumsum,maxCumsum),
#      xlab="genome index",
#      ylab="cumsum FC product window",
#      main=paste0("all DS - n = ", nDS)
# )
# mtext(side=3, text=paste0("(", minQt, "-", maxQt, "quantiles)"))
# polygon( c(1:ncol(cs_genome_qts), rev(1:ncol(cs_genome_qts))), c(cs_genome_qts[1,], rev(cs_genome_qts[2,])), col=genomeCol, 
#          density = c(20, 60), angle = c(-45, -45))
# polygon( c(1:ncol(cs_genome_tads_qts), rev(1:ncol(cs_genome_tads_qts))), c(cs_genome_tads_qts[1,], rev(cs_genome_tads_qts[2,])), col=genomeTadCol, 
#          density = c(20, 60), angle = c(-45, -45))
# foo <- dev.off()
# cat(paste0("... written: ", outFile,"\n"))

minCumSumMean <- min(c(cs_genome_mean, cs_genome_tads_mean, cs_random005_mean, cs_random095_mean))
maxCumSumMean <- max(c(cs_genome_mean, cs_genome_tads_mean, cs_random005_mean, cs_random095_mean))

outFile <- file.path(outFolder, paste0("all_ds_cumsumcurves_meanDS.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(NULL,
     #      type="l", col=genomeCol,
     xlim=c(0,minLength),
     ylim=c(minCumSumMean,maxCumSumMean),
     xlab="genome index",
     ylab="cumsum FC product window",
     main=paste0("all DS - n = ", nDS)
)
mtext(side=3, text=paste0("(all ds and ", minQt, "-", maxQt, " quantiles mean values)"))
lines(x=1:length(cs_genome_mean), y=cs_genome_mean, type="l", col=genomeCol)
lines(x=1:length(cs_genome_tads_mean), y=cs_genome_tads_mean, type="l", col=genomeTadCol)
polygon( c(1:length(cs_random005_mean), rev(1:length(cs_random095_mean))), c(cs_random005_mean, rev(cs_random095_mean)), col=randomQtCol, 
         density = c(20, 60), angle = c(-45, -45))
legend("topleft", lty=1, legend=c("no TADs", "with TADs", "with random TADs"), col = c(genomeCol, genomeTadCol, randomQtCol), bty="n", cex=1)

foo <- dev.off()
cat(paste0("... written: ", outFile,"\n"))






######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





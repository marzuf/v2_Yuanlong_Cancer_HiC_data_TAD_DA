startTime <- Sys.time()
cat(paste0("> Rscript ahlfors_FC_test.R\n"))

plotType <- "png"
myHeight  <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight


myHeightGG <- 7
myWidthGG <- myHeightGG*1.2


# Rscript ahlfors_flips_test.R

script_name <- "ahlfors_flips_test.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS=F
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

require(ggpubr)

buildTable <- TRUE

hicds = "Panc1_rep12_40kb"
exprds = "TCGApaad_wt_mutKRAS"

mainFolder <- file.path(".")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"


randomQtCol <- "lightblue"
nRandom <- 1000
set.seed(20190917)
genomeCol <- "black"
genomeTadCol <- "red"

minQt <- 0.05
maxQt <- 0.95

outFolder <- file.path("AHLFORS_FLIPS_TEST")
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
      all_sumFC <- foreach(i = 2:c(nrow(g2t_dt))) %dopar% {
        chromo_i <- g2t_dt$chromo[i]
        chromo_prev <- g2t_dt$chromo[i-1]
        
        if( chromo_i != chromo_prev) return(NULL)
        entrez_i <- g2t_dt$entrezID[i]
        stopifnot(entrez_i %in% names(gene_logFC)); stopifnot(entrez_i %in% names(gene_tad))
        entrez_prev <- g2t_dt$entrezID[i-1]
        stopifnot(entrez_prev %in% names(gene_logFC)); stopifnot(entrez_prev %in% names(gene_tad))
        
        
        fc_i <- gene_logFC[paste0(entrez_i)]
        fc_prev <- gene_logFC[paste0(entrez_prev)]
        
        
        tad_i <- gene_tad[paste0(entrez_i)]
        tad_prev <- gene_tad[paste0(entrez_prev)]
        
        
        
        i_sign_flip <- as.numeric(sign(fc_i) != sign(fc_prev))
        
        if(tad_i == tad_prev ) {
          i_tad_sign_flip <- i_sign_flip
        } else {
          i_tad_sign_flip <- NULL
        }
        list(
          i_sign_flip = i_sign_flip,
          i_tad_sign_flip = i_tad_sign_flip)
      }
      
      all_signFlip_values <- unlist(lapply(all_sumFC, function(x)x[["i_sign_flip"]]))
      cumsum_genome <- cumsum(all_signFlip_values)
      
      all_tad_signFlip_values <- unlist(lapply(all_sumFC, function(x)x[["i_tad_sign_flip"]]))
      cumsum_genome_tads <- cumsum(all_tad_signFlip_values)
      
      nBoundaries <-  length(cumsum_genome) - length(cumsum_genome_tads)
      
      # all_random_cumsum <- foreach(irand = 1:nRandom, .combine='rbind') %dopar% {
      #   random_null_idx <- sample(1:length(cumsum_genome), size=nBoundaries, replace = FALSE)
      #   random_values <- all_signFlip_values[-random_null_idx]
      #   stopifnot(length(random_values) == length(all_tad_signFlip_values))
      #   cumsum_rand <- cumsum(random_values)
      #   #lines(x=1:length(cumsum_rand), y=cumsum_rand, type="l", col="grey")
      #   cumsum_rand
      # }
      all_random <- foreach(irand = 1:nRandom) %dopar% {
        random_null_idx <- sample(1:length(cumsum_genome), size=nBoundaries, replace = FALSE)
        random_values <- all_signFlip_values[-random_null_idx]
        stopifnot(length(random_values) == length(all_tad_signFlip_values))
        cumsum_rand <- cumsum(random_values)
        #lines(x=1:length(cumsum_rand), y=cumsum_rand, type="l", col="grey")
        list(cumsum_rand=cumsum_rand,random_values=random_values)
      }
      
      all_random_cumsum <- do.call(rbind, lapply(all_random, function(x)x[["cumsum_rand"]]))
      all_random_values <- do.call(rbind, lapply(all_random, function(x)x[["random_values"]]))
      
      
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
        all_signFlip_values=all_signFlip_values,
        all_tad_signFlip_values=all_tad_signFlip_values,
        all_random_values = all_random_values,
        cumsum_genome=cumsum_genome,
        cumsum_genome_tads=cumsum_genome_tads,
        cumsum_random = all_random_cumsum,
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
  # outFile= "AHLFORS_FLIPS_TEST/all_cumsum_data.Rdata"
  outFile <- file.path(outFolder, paste0("all_cumsum_data.Rdata"))
  cat(paste0("... load data\n"))
  all_cumsum_data <- get(load(outFile))
  
}

all_cumsum_data_ul <- unlist(all_cumsum_data, recursive=FALSE)

all_nFlips_genome <- unlist(lapply(lapply(all_cumsum_data_ul,function(x) x[["all_signFlip_values"]]), sum))
all_nFlips_genome_tad <- unlist(lapply(lapply(all_cumsum_data_ul,function(x) x[["all_tad_signFlip_values"]]), sum))
all_nFlips_random <- unlist(lapply(lapply(all_cumsum_data_ul, function(x) x[["all_random_values"]]), function(dt)apply(dt,1,sum)))


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
     ylab="cumsum sign flip product window",
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


all_nFlips_genome <- unlist(lapply(lapply(all_cumsum_data_ul,function(x) x[["all_signFlip_values"]]), sum))
all_nFlips_genome_tad <- unlist(lapply(lapply(all_cumsum_data_ul,function(x) x[["all_tad_signFlip_values"]]), sum))
all_nFlips_random <- unlist(lapply(lapply(all_cumsum_data_ul, function(x) x[["all_random_values"]]), function(dt)apply(dt,1,sum)))

plot_dt <- rbind(
  rbind(
  data.frame(
    value_type = "genome",
    flip_values = c(all_nFlips_genome),
    stringsAsFactors = FALSE
  ),
  data.frame(
    value_type = "with_tads",
    flip_values = c(all_nFlips_genome_tad),
    stringsAsFactors = FALSE
  )
),
data.frame(
  value_type = "random_tads",
  flip_values = c(all_nFlips_random),
  stringsAsFactors = FALSE
)
)

plot_dt <- plot_dt[!is.null(plot_dt$flip_values),]

plot_dt$value_type <- factor(plot_dt$value_type, levels=c("genome", "with_tads", "random_tads"))


p <- ggdensity(plot_dt, 
               x = "flip_values", 
            xlab="sum # flips",
               color = "value_type", fill = "value_type",
               add = "mean", rug = TRUE,
               palette = c(genomeCol, genomeTadCol, randomQtCol))

outFile <- file.path(outFolder, paste0("all_ds_values_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))



all_meanFlips_genome <- unlist(lapply(lapply(all_cumsum_data_ul,function(x) x[["all_signFlip_values"]]), mean))
all_meanFlips_genome_tad <- unlist(lapply(lapply(all_cumsum_data_ul,function(x) x[["all_tad_signFlip_values"]]), mean))
all_meanFlips_random <- unlist(lapply(lapply(all_cumsum_data_ul, function(x) x[["all_random_values"]]), function(dt)apply(dt,1,mean)))

plot_dt <- rbind(
  rbind(
  data.frame(
    value_type = "genome",
    flip_values = c(all_meanFlips_genome),
    stringsAsFactors = FALSE
  ),
  data.frame(
    value_type = "with_tads",
    flip_values = c(all_meanFlips_genome_tad),
    stringsAsFactors = FALSE
  )
),
data.frame(
  value_type = "random_tads",
  flip_values = c(all_meanFlips_random),
  stringsAsFactors = FALSE
)
)

plot_dt <- plot_dt[!is.null(plot_dt$flip_values),]

plot_dt$value_type <- factor(plot_dt$value_type, levels=c("genome", "with_tads", "random_tads"))


p <- ggdensity(plot_dt, 
               x = "flip_values", 
            xlab="mean # flips",
               color = "value_type", fill = "value_type",
               add = "mean", rug = TRUE,
               palette = c(genomeCol, genomeTadCol, randomQtCol))

outFile <- file.path(outFolder, paste0("all_ds_mean_values_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))

######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))






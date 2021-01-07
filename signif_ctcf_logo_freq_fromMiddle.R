
# Rscript signif_ctcf_logo_freq_fromMiddle.R

set.seed(21122020)

outFolder <- file.path("SIGNIF_CTCF_LOGO_FREQ_FROMMIDDLE")
dir.create(outFolder, recursive = TRUE)

# source("http://bioconductor.org/biocLite.R")
# BiocManager::install("Logolas")
library("Logolas")
library("Biostrings")
library("msa")
library("ggseqlogo")
require(ggplot2)
library(grid)
library(gridExtra)
library(stringr)
library(foreach)

plotType <- "svg"
myHeightGG <- 5
myWidthGG <- 6

# dnaSet = DNAStringSet(c("AACCTT","CCGGTTTT","AAAGGGTTT"))
# res = msa(dnaSet,method="ClustalOmega")
# print(res)
# 
# myctcf_seqs <- c(">>><<",">>>><","><><>>")
# myctcf_seqs_hack <- gsub("<", "T", gsub(">", "A", myctcf_seqs))
# dnaSet = DNAStringSet(myctcf_seqs_hack)
# res = msa(dnaSet,method="ClustalOmega")
# print(res)
# aligned_ctcf_sequences <- as.character(res)
# logomaker(aligned_ctcf_sequences, type = "Logo")
# aligned_ctcf_sequences_hack <- gsub("T", "<", gsub("A", ">", aligned_ctcf_sequences))
# logomaker(aligned_ctcf_sequences_hack, type = "Logo")

logFile <- file.path(outFolder, "signif_ctcf_logo_freq_fromMiddle_logFile.txt")
# system(paste0("rm -f ))
sink(logFile)

ctcf2tad_dt <- get(load("CTCF_AND_DA_ALLDS/ctcf2tad_dt.Rdata"))
all_result_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

pthresh <- 0.01

signif_dt <- all_result_dt[all_result_dt$adjPvalComb <= pthresh,]
nonsignif_dt <- all_result_dt[all_result_dt$adjPvalComb > pthresh,]

nSignif <- sum(all_result_dt$adjPvalComb <= pthresh)
random_dt <- all_result_dt[sample(1:nrow(all_result_dt), size=nSignif),]
stopifnot(!is.na(random_dt))
stopifnot(nrow(random_dt) == nrow(signif_dt))

# prepare random same size dist
############################################################################################################
# random v2 -> keep size distribution ! 
all_ctcf2tad_dt <- merge(all_result_dt,ctcf2tad_dt, by=c("hicds", "region"), all=FALSE , suffixes=c("_TAD", "_CTCF"))
all_ctcf2tad_dt$chromo <- gsub("(chr.+)_.+", "\\1", all_ctcf2tad_dt$region)
stopifnot(all_ctcf2tad_dt$chromo %in% paste0("chr", 1:22))
all_ctcf2tad_dt <- all_ctcf2tad_dt[order(all_ctcf2tad_dt$chromo, all_ctcf2tad_dt$start_CTCF, all_ctcf2tad_dt$end_CTCF),]

all_ctcf2tad_dt$ctcf_midPos <- 0.5*(all_ctcf2tad_dt$start_CTCF + all_ctcf2tad_dt$end_CTCF)
all_ctcf2tad_dt$tad_midPos <- 0.5*(all_ctcf2tad_dt$start_TAD + all_ctcf2tad_dt$end_TAD)

left_ctcf2tad_dt <- all_ctcf2tad_dt[all_ctcf2tad_dt$ctcf_midPos <= all_ctcf2tad_dt$tad_midPos,]
right_ctcf2tad_dt <- all_ctcf2tad_dt[all_ctcf2tad_dt$ctcf_midPos > all_ctcf2tad_dt$tad_midPos,]
stopifnot(nrow(right_ctcf2tad_dt) + nrow(left_ctcf2tad_dt) == nrow(all_ctcf2tad_dt))

rm("all_ctcf2tad_dt")

all_dirs <- c("left", "right")
tad_dir="left"

for(tad_dir in all_dirs) {
  
  
  motif_agg_dt <- aggregate(orientation~hicds + exprds + region + adjPvalComb, data=eval(parse(text = paste0(tad_dir, "_ctcf2tad_dt"))), 
                            FUN=function(x) paste0(rev(x), collapse=""))  #### here use "rev" for fromMiddle version
  motif_agg_dt$regionID <- file.path(motif_agg_dt$hicds, motif_agg_dt$exprds, motif_agg_dt$region)
  stopifnot(!duplicated(motif_agg_dt$regionID))
  # nrow(motif_agg_dt)
  
  signif_motif_agg_dt <- motif_agg_dt[motif_agg_dt$adjPvalComb <= pthresh,]
  signif_motif_freqDT <- data.frame(table(signif_motif_agg_dt$orientation))
  colnames(signif_motif_freqDT)[1] <- "Motif"
  signif_motif_freqDT <- signif_motif_freqDT[order(signif_motif_freqDT$Freq, decreasing=T),]
  signif_motif_freqDT$FreqRatio <- round(signif_motif_freqDT$Freq/sum(signif_motif_freqDT$Freq),4)
  cat(paste0(tad_dir, " - ", "head(signif_motif_freqDT)", "\n"))
  # head(signif_motif_freqDT)
  write.table(signif_motif_freqDT[1:10,], quote=F, sep="\t", row.names=F)
  
  nonsignif_motif_agg_dt <- motif_agg_dt[motif_agg_dt$adjPvalComb > pthresh,]
  nonsignif_motif_freqDT <- data.frame(table(nonsignif_motif_agg_dt$orientation))
  colnames(nonsignif_motif_freqDT)[1] <- "Motif"
  nonsignif_motif_freqDT <- nonsignif_motif_freqDT[order(nonsignif_motif_freqDT$Freq, decreasing=T),]
  nonsignif_motif_freqDT$FreqRatio <- round(nonsignif_motif_freqDT$Freq/sum(nonsignif_motif_freqDT$Freq),4)
  cat(paste0(tad_dir, " - ", "head(nonsignif_motif_freqDT)", "\n"))
  # head(nonsignif_motif_freqDT)
  write.table(nonsignif_motif_freqDT[1:10,], quote=F, sep="\t", row.names=F)
  
  
  
  seqmotif_agg_dt <- aggregate(orientation~hicds + exprds + region + adjPvalComb, data=eval(parse(text = paste0(tad_dir, "_ctcf2tad_dt"))), 
                               FUN=function(x)  paste0(rle(x)$values, collapse=""))
  seqmotif_agg_dt$regionID <- file.path(seqmotif_agg_dt$hicds, seqmotif_agg_dt$exprds, seqmotif_agg_dt$region)
  stopifnot(!duplicated(seqmotif_agg_dt$regionID))
  # nrow(seqmotif_agg_dt)
  
  
  signif_seqmotif_agg_dt <- seqmotif_agg_dt[seqmotif_agg_dt$adjPvalComb <= pthresh,]
  signif_seqmotif_freqDT <- data.frame(table(signif_seqmotif_agg_dt$orientation))
  colnames(signif_seqmotif_freqDT)[1] <- "Motif"
  signif_seqmotif_freqDT <- signif_seqmotif_freqDT[order(signif_seqmotif_freqDT$Freq, decreasing=T),]
  signif_seqmotif_freqDT$FreqRatio <- round(signif_seqmotif_freqDT$Freq/sum(signif_seqmotif_freqDT$Freq),4)
  cat(paste0(tad_dir, " - ", "head(signif_seqmotif_freqDT)", "\n"))
  # head(signif_seqmotif_freqDT)
  write.table(signif_seqmotif_freqDT[1:10,], quote=F, sep="\t", row.names=F)
  
  nonsignif_seqmotif_agg_dt <- seqmotif_agg_dt[seqmotif_agg_dt$adjPvalComb > pthresh,]
  nonsignif_seqmotif_freqDT <- data.frame(table(nonsignif_seqmotif_agg_dt$orientation))
  colnames(nonsignif_seqmotif_freqDT)[1] <- "Motif"
  nonsignif_seqmotif_freqDT <- nonsignif_seqmotif_freqDT[order(nonsignif_seqmotif_freqDT$Freq, decreasing=T),]
  nonsignif_seqmotif_freqDT$FreqRatio <- round(nonsignif_seqmotif_freqDT$Freq/sum(nonsignif_seqmotif_freqDT$Freq),4)
  cat(paste0(tad_dir, " - ", "head(nonsignif_seqmotif_freqDT)", "\n"))
  # head(nonsignif_seqmotif_freqDT)
  write.table(nonsignif_seqmotif_freqDT[1:10,], quote=F, sep="\t", row.names=F)
  
  all_dts <- c("signif_seqmotif_agg_dt", "nonsignif_seqmotif_agg_dt", "signif_motif_agg_dt", "nonsignif_motif_agg_dt")
  
  # iterate for each # of motifs
  # do the MSA 
  # 
                    # in_dt <- signif_motif_agg_dt
                    # dt=all_dts[1]
                    # for(dt in all_dts) {
                    #   
                    #   in_dt <- eval(parse(text = dt))
                    #   
                    #   in_dt$nSites <- str_count(in_dt$orientation)
                    #   
                    #   all_nsites <- unique(in_dt$nSites)
                    #   nsites  = all_nsites[1]
                    #   sub_dt=in_dt[in_dt$nSites==nsites,]
                    #   
                    #   
                    #   out_dt <- foreach(nsites = all_nsites, .combine='rbind') %do% {
                    #     
                    #     sub_dt <- in_dt[in_dt$nSites==nsites,,drop=FALSE]
                    #     if(nrow(sub_dt) == 1 | nsites==1) {
                    #       aligned_ctcf_seqs <- sub_dt$orientation
                    #       
                    #     } else {
                    #       ctcf_seqs <- sub_dt$orientation
                    #       ctcf_seqs_hack <- gsub("<", "T", gsub(">", "A", ctcf_seqs))
                    #       dnaSet <- DNAStringSet(ctcf_seqs_hack)
                    #       res_aligned_seqs <- msa(dnaSet,method="ClustalOmega")
                    #       print(res_aligned_seqs)
                    #       aligned_ctcf_seqs_hack <- as.character(res_aligned_seqs)
                    #       # logomaker(aligned_ctcf_sequences, type = "Logo")
                    #       aligned_ctcf_seqs <- gsub("T", "<", gsub("A", ">", aligned_ctcf_seqs_hack))
                    #       
                    #     }
                    #     
                    #     msa_freqDT <- data.frame(table(aligned_ctcf_seqs))
                    #     colnames(msa_freqDT)[1] <- "Motif"
                    #     msa_freqDT <- msa_freqDT[order(msa_freqDT$Freq, decreasing=T),]
                    #     msa_freqDT$FreqRatio <- msa_freqDT$Freq/sum(msa_freqDT$Freq)
                    #     # head(msa_freqDT)
                    #     msa_freqDT$nSites <- nsites
                    #     
                    #     cat(paste0(tad_dir, " - ", dt, " - ", nsites, " - ", "head(msa_freqDT)", "\n"))
                    #     # head(msa_freqDT)
                    #     write.table(msa_freqDT[1:10,], quote=F, sep="\t", row.names=F)
                    #     
                    #     msa_freqDT
                    #   }
                    #   out_dt <- out_dt[order(out_dt$Freq, out_dt$FreqRatio, decreasing=T),]
                    # }
  
  in_dt <- signif_motif_agg_dt
  dt=all_dts[1]
  for(dt in all_dts) {
    
    in_dt <- eval(parse(text = dt))
    
    in_dt$nSites <- str_count(in_dt$orientation)
    
    all_nsites <- unique(in_dt$nSites)
    nsites  = all_nsites[1]
    sub_dt=in_dt[in_dt$nSites==nsites,]
    
    
    out_dt <- foreach(nsites = all_nsites, .combine='rbind') %do% {
      
      sub_dt <- in_dt[in_dt$nSites==nsites,,drop=FALSE]

      
      sub_freqDT <- data.frame(table(sub_dt$orientation))
      colnames(sub_freqDT)[1] <- "Motif"
      sub_freqDT <- sub_freqDT[order(sub_freqDT$Freq, decreasing=T),]
      sub_freqDT$FreqRatio <- round(sub_freqDT$Freq/sum(sub_freqDT$Freq),4)
      cat(paste0(tad_dir, " - ", dt, " - ", nsites, " - ", "head(sub_freqDT)", "\n"))      # head(signif_seqmotif_freqDT)
      write.table(sub_freqDT[1:10,], quote=F, sep="\t", row.names=F)
      
      sub_freqDT
    }
    out_dt <- out_dt[order(out_dt$Freq, out_dt$FreqRatio, decreasing=T),]
  }
  
}





sink()

cat(paste0("... written: ", logFile, "\n"))

stop("-ok\n")













all_agg_dt <- aggregate(orientation~hicds + exprds + region + adjPvalComb, data=all_ctcf2tad_dt, FUN=length)
all_agg_dt$regionID <- file.path(all_agg_dt$hicds, all_agg_dt$exprds, all_agg_dt$region)
nrow(all_agg_dt)
colnames(all_agg_dt)[colnames(all_agg_dt) == "orientation"] <- "nBS"

signif_agg_dt <- all_agg_dt[all_agg_dt$adjPvalComb <= pthresh,]
nonsignif_agg_dt <- all_agg_dt[all_agg_dt$adjPvalComb > pthresh,]





all_sizes_dt <- aggregate(region~nBS, data=signif_agg_dt, FUN=length)
stopifnot(!duplicated(all_sizes_dt$nBS))

sampledTADs <- foreach(i = 1:nrow(all_sizes_dt), .combine='c') %dopar% {
  
  nBS <- all_sizes_dt$nBS[i]
  nToSample <- all_sizes_dt$region[i]  
  
  bag_dt <- nonsignif_agg_dt[nonsignif_agg_dt$nBS == nBS,]
  stopifnot(nrow(bag_dt) > nToSample)
  
  bag_dt$regionID[sample(c(1:nrow(bag_dt)), size = nToSample)]
}
random_same_size_dt <- all_result_dt
random_same_size_dt$region_id <- file.path(random_same_size_dt$hicds, random_same_size_dt$exprds, random_same_size_dt$region)

stopifnot(length(sampledTADs) == nrow(signif_agg_dt))
stopifnot(!duplicated(sampledTADs) )
stopifnot(sampledTADs %in% file.path(nonsignif_dt$hicds, nonsignif_dt$exprds, nonsignif_dt$region))
stopifnot(sampledTADs %in% random_same_size_dt$region_id)
random_same_size_dt <- random_same_size_dt[random_same_size_dt$region_id %in% sampledTADs,]
random_same_size_dt$region_id <- NULL

stopifnot(nrow(random_same_size_dt) == nrow(signif_agg_dt))


for(dataType in c("signif.", "non signif.", "random", "random vSameSizeDist")) {
  # for(dataType in c("random")) {
  
  cat(paste0("... START\t", dataType, "\n"))
  
  
  # > summary(signif_agg_dt$nBS)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 1.00    6.00    9.00   10.32   14.00   40.00 
  # 
  
  if(dataType == "signif.") {
    input_dt <- signif_dt  
    # thresh1 <- 6 # 8
    # thresh2 <- 14 # 12
    thresh1 <-  8
    thresh2 <-  12
    dataTypeOut <- "signif"
  } else if(dataType == "non signif."){
    input_dt <- nonsignif_dt  
    # thresh1 <- 6 # 8
    # thresh2 <- 14 # 12
    thresh1 <-  8
    thresh2 <- 12
    dataTypeOut <- "nonsignif"
  } else if(dataType == "random"){
    input_dt <- random_dt  
    # thresh1 <- 6 # 8
    # thresh2 <- 14 # 12
    thresh1 <-  8
    thresh2 <-  12
    dataTypeOut <- "random"
  }  else if(dataType == "random vSameSizeDist"){
    input_dt <- random_same_size_dt  
    # thresh1 <- 6 # 8
    # thresh2 <- 14 # 12
    thresh1 <-  8
    thresh2 <-  12
    dataTypeOut <- "randomVsameSizeDist"
  }else{
    stop("error\n")
  }
  
  plotTit <- paste0(dataType, " TADs - n=", nrow(input_dt))
  
  input_ctcf2tad_dt <- merge(input_dt,ctcf2tad_dt, by=c("hicds", "region"), all=FALSE , suffixes=c("_TAD", "_CTCF"))
  input_ctcf2tad_dt$chromo <- gsub("(chr.+)_.+", "\\1", input_ctcf2tad_dt$region)
  stopifnot(input_ctcf2tad_dt$chromo %in% paste0("chr", 1:22))
  input_ctcf2tad_dt <- input_ctcf2tad_dt[order(input_ctcf2tad_dt$chromo, input_ctcf2tad_dt$start_CTCF, input_ctcf2tad_dt$end_CTCF),]
  
  agg_dt <- aggregate(orientation~hicds + exprds + region, data=input_ctcf2tad_dt, FUN=function(x) paste0(x, collapse=""))
  agg_dt$regionID <- file.path(agg_dt$hicds, agg_dt$exprds, agg_dt$region)
  stopifnot(!duplicated(agg_dt$regionID))
  nrow(agg_dt)
  # 1442
  # NS: 99203
  
  if(dataType == "signif.") {
    tmp <- agg_dt
    tmp$nMotifs <- nchar(tmp$orientation)
    tmp <- tmp[order(-tmp$nMotifs, tmp$hicds, tmp$exprds),]
    outFile <- file.path(outFolder, paste0("CTCForientationSeq_", dataTypeOut, "_dt.txt"))
    write.table(tmp[,c("regionID", "orientation", "nMotifs")], file=outFile, col.names = F, row.names=F)
    cat(paste0("... written: ", outFile, "\n"))
    stop("-ok")
  }
  
  length_agg_dt <- aggregate(orientation~hicds + exprds + region, data=input_ctcf2tad_dt, FUN=length)
  length_agg_dt$regionID <- file.path(length_agg_dt$hicds, length_agg_dt$exprds, length_agg_dt$region)
  
  outFile <- file.path(outFolder, paste0("nbrCTCFbindingsites_", dataTypeOut, "TADs_density.", plotType))
  do.call(plotType, list(outFile, height=myHeightGG, width=myWidthGG))
  plot(density(length_agg_dt$orientation), main=paste0("# CTCF BS by TAD - ", dataType))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  ########################### take only those between 8 and 12 for the signif
  keepTADs <- length_agg_dt$regionID[length_agg_dt$orientation >= thresh1 & length_agg_dt$orientation <= thresh2]
  length(keepTADs)
  # 474
  
  subTit <- paste0("# CTCF BS >=", thresh1, " & <=", thresh2, " (n=", length(keepTADs), ")")
  
  sub_agg_dt <- agg_dt[agg_dt$regionID %in% keepTADs,]
  stopifnot(nrow(sub_agg_dt) == length(keepTADs))
  
  
  myctcf_seqs <- sub_agg_dt$orientation
  myctcf_seqs_hack <- gsub("<", "T", gsub(">", "A", myctcf_seqs))
  dnaSet = DNAStringSet(myctcf_seqs_hack)
  res = msa(dnaSet,method="ClustalOmega")
  print(res)
  aligned_ctcf_sequences <- as.character(res)
  # logomaker(aligned_ctcf_sequences, type = "Logo")
  aligned_ctcf_sequences_hack <- gsub("T", "<", gsub("A", ">", aligned_ctcf_sequences))
  
  outFile <- file.path(outFolder, paste0(dataTypeOut, "TADs_alignedSeqLogo_nSites_", thresh1, "_", thresh2, "_logo.", plotType))
  do.call(plotType, list(outFile, height=myHeightGG, width=myWidthGG))
  logomaker(aligned_ctcf_sequences_hack, type = "Logo")
  grid.text( plotTit, x = unit(0.5, "npc"), y = unit(0.9, "npc"),gp = gpar(fontface = "bold"))
  grid.text( subTit, x = unit(0.5, "npc"), y = unit(0.8, "npc"),gp = gpar(fontface = "italic"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  aligned_ctcf_sequences_hack2 <- gsub("T", "R", gsub("A", "F", aligned_ctcf_sequences))
  outFile <- file.path(outFolder, paste0(dataTypeOut, "TADs_alignedSeqLogo_nSites_", thresh1, "_", thresh2, "_logo_hack2.", plotType))
  do.call(plotType, list(outFile, height=myHeightGG, width=myWidthGG))
  logomaker(aligned_ctcf_sequences_hack2, type = "Logo")
  grid.text( plotTit, x = unit(0.5, "npc"), y = unit(0.9, "npc"),gp = gpar(fontface = "bold"))
  grid.text( subTit, x = unit(0.5, "npc"), y = unit(0.8, "npc"),gp = gpar(fontface = "italic"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  aligned_ctcf_sequences_hack2 <- gsub("T", "R", gsub("A", "F", aligned_ctcf_sequences))
  p <- ggseqlogo(aligned_ctcf_sequences_hack2, seq_type="other", namespace=c("R", "F")) + 
    ggtitle(plotTit, subtitle=subTit) + theme(plot.title = element_text(face="bold", hjust=0.5),
                                              plot.subtitle = element_text(face="bold", hjust=0.5)
                                              )
  outFile <- file.path(outFolder, paste0(dataTypeOut, "TADs_alignedSeqLogo_nSites_", thresh1, "_", thresh2, "_ggseq_hack2.", plotType))
  ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  ########################### take only smaller than thresh1
  
  
  
  keepTADs <- length_agg_dt$regionID[length_agg_dt$orientation < thresh1 ]
  length(keepTADs)
  # 547
  
  subTit <- paste0("# CTCF BS <", thresh1, " (n=", length(keepTADs), ")")

  sub_agg_dt <- agg_dt[agg_dt$regionID %in% keepTADs,]
  stopifnot(nrow(sub_agg_dt) == length(keepTADs))

  myctcf_seqs <- sub_agg_dt$orientation
  myctcf_seqs_hack <- gsub("<", "T", gsub(">", "A", myctcf_seqs))
  dnaSet = DNAStringSet(myctcf_seqs_hack)
  res = msa(dnaSet,method="ClustalOmega")
  print(res)
  aligned_ctcf_sequences <- as.character(res)
  logomaker(aligned_ctcf_sequences, type = "Logo")

  aligned_ctcf_sequences_hack <- gsub("T", "<", gsub("A", ">", aligned_ctcf_sequences))

  outFile <- file.path(outFolder, paste0(dataTypeOut, "TADs_alignedSeqLogo_nSites_smaller", thresh1, "_logo.", plotType))
  do.call(plotType, list(outFile, height=myHeightGG, width=myWidthGG))
  logomaker(aligned_ctcf_sequences_hack, type = "Logo")
  grid.text( plotTit, x = unit(0.5, "npc"), y = unit(0.9, "npc"),gp = gpar(fontface = "bold"))
  grid.text( subTit, x = unit(0.5, "npc"), y = unit(0.8, "npc"),gp = gpar(fontface = "italic"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  aligned_ctcf_sequences_hack2 <- gsub("T", "R", gsub("A", "F", aligned_ctcf_sequences))
  outFile <- file.path(outFolder, paste0(dataTypeOut, "TADs_alignedSeqLogo_nSites_smaller", thresh1, "_logo_hack2.", plotType))
  do.call(plotType, list(outFile, height=myHeightGG, width=myWidthGG))
  logomaker(aligned_ctcf_sequences_hack2, type = "Logo")
  grid.text( plotTit, x = unit(0.5, "npc"), y = unit(0.9, "npc"),gp = gpar(fontface = "bold"))
  grid.text( subTit, x = unit(0.5, "npc"), y = unit(0.8, "npc"),gp = gpar(fontface = "italic"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  aligned_ctcf_sequences_hack2 <- gsub("T", "R", gsub("A", "F", aligned_ctcf_sequences))
  p <- ggseqlogo(aligned_ctcf_sequences_hack2, seq_type="other", namespace=c("R", "F")) + 
    ggtitle(plotTit, subtitle=subTit) + theme(plot.title = element_text(face="bold", hjust=0.5),
                                              plot.subtitle = element_text(face="bold", hjust=0.5)
    )
  outFile <- file.path(outFolder, paste0(dataTypeOut, "TADs_alignedSeqLogo_nSites_smaller", thresh1, "_ggseq_hack2.", plotType))
  ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  

  ## !! position 3 and 11 !!!

  ########################### take only bigger than thresh2

  keepTADs <- length_agg_dt$regionID[length_agg_dt$orientation > thresh2]
  length(keepTADs)
  # 421
  
  subTit <- paste0("# CTCF BS >", thresh2, " (n=", length(keepTADs), ")")

  sub_agg_dt <- agg_dt[agg_dt$regionID %in% keepTADs,]
  stopifnot(nrow(sub_agg_dt) == length(keepTADs))


  myctcf_seqs <- sub_agg_dt$orientation
  myctcf_seqs_hack <- gsub("<", "T", gsub(">", "A", myctcf_seqs))
  dnaSet = DNAStringSet(myctcf_seqs_hack)
  res = msa(dnaSet,method="ClustalOmega")
  print(res)
  aligned_ctcf_sequences <- as.character(res)
  logomaker(aligned_ctcf_sequences, type = "Logo")

  aligned_ctcf_sequences_hack <- gsub("T", "<", gsub("A", ">", aligned_ctcf_sequences))

  outFile <- file.path(outFolder, paste0(dataTypeOut, "TADs_alignedSeqLogo_nSites_bigger", thresh2, "_logo.", plotType))
  do.call(plotType, list(outFile, height=myHeightGG, width=myWidthGG))
  logomaker(aligned_ctcf_sequences_hack, type = "Logo")
  grid.text( plotTit, x = unit(0.5, "npc"), y = unit(0.9, "npc"),gp = gpar(fontface = "bold"))
  grid.text( subTit, x = unit(0.5, "npc"), y = unit(0.8, "npc"),gp = gpar(fontface = "italic"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  aligned_ctcf_sequences_hack2 <- gsub("T", "R", gsub("A", "F", aligned_ctcf_sequences))
  outFile <- file.path(outFolder, paste0(dataTypeOut, "TADs_alignedSeqLogo_nSites_bigger", thresh2, "_logo_hack2.", plotType))
  do.call(plotType, list(outFile, height=myHeightGG, width=myWidthGG))
  logomaker(aligned_ctcf_sequences_hack2, type = "Logo")
  grid.text( plotTit, x = unit(0.5, "npc"), y = unit(0.9, "npc"),gp = gpar(fontface = "bold"))
  grid.text( subTit, x = unit(0.5, "npc"), y = unit(0.8, "npc"),gp = gpar(fontface = "italic"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  aligned_ctcf_sequences_hack2 <- gsub("T", "R", gsub("A", "F", aligned_ctcf_sequences))
  p <- ggseqlogo(aligned_ctcf_sequences_hack2, seq_type="other", namespace=c("R", "F")) + 
    ggtitle(plotTit, subtitle=subTit) + theme(plot.title = element_text(face="bold", hjust=0.5),
                                              plot.subtitle = element_text(face="bold", hjust=0.5)
    )
  outFile <- file.path(outFolder, paste0(dataTypeOut, "TADs_alignedSeqLogo_nSites_bigger", thresh2, "_ggseq_hack2.", plotType))
  ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
}


# Rscript cmp_synBlocks.R

infercars_dt <- read.delim("inferCARs_data/Orthology.Blocks_processed_hg19.txt", stringsAsFactors = FALSE, header=TRUE)
procars_dt <- read.delim("procars_orthology_blocks_processsed.txt", stringsAsFactors = FALSE, header=TRUE)

infercars_dt <- infercars_dt[infercars_dt$genome == "hg19",]
procars_dt <- procars_dt[procars_dt$genome == "homo_sapiens",]

stopifnot(!duplicated(infercars_dt$blockID))
stopifnot(!duplicated(procars_dt$blockID))
# nrow(infercars_dt)
# [1] 3135
# > nrow(procars_dt)
# [1] 689

require(GenomicRanges)
require(ggsci)

plotType <- "png"
myWidth <- 600
myHeight <- 400
plotCex <- 1.4

infer_col <- pal_d3()(2)[1]
pro_col <- pal_d3()(2)[2]

outFolder <- file.path("CMP_SYNBLOCKS")
dir.create(outFolder, recursive = TRUE)

logFile <- file.path(outFolder, "cmp_syn_blocks.txt")
file.remove(logFile)
file.create(logFile)

infercars_GR <- GRanges(seqnames = infercars_dt$chromo, ranges = IRanges(start=infercars_dt$start, 
                                                                     end=infercars_dt$end, 
                                                                     names=infercars_dt$blockID))

procars_GR <- GRanges(seqnames = procars_dt$chromo, ranges = IRanges(start=procars_dt$start, 
                                                                         end=procars_dt$end, 
                                                                         names=procars_dt$blockID))
overlap_GR <- findOverlaps(query=infercars_GR, subject=procars_GR)

queryID <- names(infercars_GR[queryHits(overlap_GR)])
refID <- names(procars_GR[subjectHits(overlap_GR)])

cars_overlap_dt <- data.frame(
  queryID_inferCARs = queryID,
  refID_proCARs = refID,
  width_inferCARs = width(infercars_GR[queryHits(overlap_GR)]),
  width_proCARs = width(procars_GR[subjectHits(overlap_GR)]),
    overlapBp = width(pintersect(procars_GR[refID], infercars_GR[queryID])),
  stringsAsFactors = FALSE)

stopifnot(cars_overlap_dt$queryID_inferCARs %in% infercars_dt$blockID)
stopifnot(cars_overlap_dt$refID_proCARs %in% procars_dt$blockID)

# inferCARs with no matching
txt <- paste0("# inferCARs with no matching:\t", sum(!infercars_dt$blockID %in% cars_overlap_dt$queryID_inferCARs), 
              "/", length(unique(infercars_dt$blockID)), "\n")
cat(paste0(txt))   
cat(txt, file = logFile, append=TRUE)

txt <- paste0("% inferCARs with no matching:\t", round(100*sum(!infercars_dt$blockID %in% cars_overlap_dt$queryID_inferCARs)/
                                                         length(unique(infercars_dt$blockID)), 2), "%\n")
cat(paste0(txt))                     
cat(txt, file = logFile, append=TRUE)

infercars_dt$width <- infercars_dt$end - infercars_dt$start
txt <- paste0("% bp inferCARs with matching:\t", round(100*sum(cars_overlap_dt$overlapBp)/sum(as.numeric(infercars_dt$width)), 2), "%\n")
cat(paste0(txt))                     
cat(txt, file = logFile, append=TRUE)



# proCARS with no matching
txt <- paste0("# proCARs with no matching:\t", sum(!procars_dt$blockID %in% cars_overlap_dt$refID_proCARs), 
              "/", length(unique(procars_dt$blockID)), "\n")
cat(paste0(txt))                                                     
cat(txt, file = logFile, append=TRUE)

txt <- paste0("% proCARs with no matching:\t", round(100*sum(!procars_dt$blockID %in% cars_overlap_dt$refID_proCARs)/
                                                         length(unique(procars_dt$blockID)), 2), "%\n")
cat(paste0(txt))                                                     
cat(txt, file = logFile, append=TRUE)

procars_dt$width <- procars_dt$end - procars_dt$start
txt <- paste0("% bp proCARs with matching:\t", round(100*sum(cars_overlap_dt$overlapBp)/sum(as.numeric(procars_dt$width)), 2), "%\n")
cat(paste0(txt))                     
cat(txt, file = logFile, append=TRUE)

# visulaization$
infercars_dt$synType <- "inferCARs"
procars_dt$synType <- "proCARs"

all_dt <- rbind(infercars_dt, procars_dt)

all_dt <- all_dt[order(all_dt$chromo, all_dt$start, all_dt$end),]

all_dt$yoffset <- ifelse(all_dt$synType == "inferCARs", -0.05, +0.05)

all_dt$segCol <- ifelse(all_dt$synType == "inferCARs", infer_col, pro_col)

all_chromo <- paste0("chr", 1:22)


outFile <- file.path(outFolder, paste0("cmp_inferCARs_proCARs_synteny_blocks.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="l")
plot(NULL,
     ylim = c(0, length(all_chromo) / 10 * 2+ 0.5), 
     xlim=c(min(all_dt$start), max(all_dt$end)),
     axes = FALSE,
     xlab="",
     ylab="",
     main = paste0("proCARs and inferCARs synteny blocks"),
     cex.main = plotCex,
     cex.axis = plotCex,
     cex.lab = plotCex
     )
axis(1)
axis(2, at = 1:length(all_chromo)/10 * 2, labels = all_chromo, las=1, lwd=0, lwd.ticks = 2)

for(i_chrom in 1:length(all_chromo)) {
  
  chr_pos <- i_chrom /10 * 2
  
  chr_dt <- all_dt[all_dt$chromo == all_chromo[i_chrom],]
  
  segments(
    x0 = chr_dt$start,
    y0 = chr_pos + chr_dt$yoffset,
    x1 = chr_dt$end,
    y1 = chr_pos + chr_dt$yoffset,
    col = chr_dt$segCol,
    lwd=2
  )
}

legend("topright",
       legend = c("proCARs", "inferCARs"),
       col = c(pro_col, infer_col),
       lty = 1,
       bty="n")


foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

cat(paste0("... written: ", logFile, "\n"))


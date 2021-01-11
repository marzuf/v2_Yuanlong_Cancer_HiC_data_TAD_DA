# CTCF_AND_FAMILY_CLUSTERS_RANDOMV2
# ctcf10000_family130000
# ctcf10000_family260000
# ctcf20000_family260000
# ctcf4000_family260000

# Rscript ctcf_clust_dist.R
require(ggpubr)

d_clustSize <- 20000

outFolder <- "CTCF_CLUST_DIST"
dir.create(outFolder, recursive = TRUE)

infile <- file.path("CTCF_AND_FAMILY_CLUSTERS_RANDOMV2", paste0("ctcf", d_clustSize, "_family260000"),
                              paste0("ctcf_bs_nanni2020_cluster", d_clustSize, "bp.bed"))

clust_dt <- read.delim(infile, header=FALSE, col.names=c("chromo", "start", "end", "orientation", "cluster"),
                       stringsAsFactors = FALSE)

clust_dt$cluster <- as.character(clust_dt$cluster)

check_dt <- aggregate(chromo~cluster, data=clust_dt, FUN=function(x) stopifnot(length(unique(x)) == 1))

aggSize_clust_dt <- aggregate(orientation~cluster, data=clust_dt, FUN=length)
colnames(aggSize_clust_dt)[2] <- "clustSize"

# x <- hist(aggSize_clust_dt$clustSize, breaks=1:10)
# x$density
# stopifnot(all(diff(x$breaks) == 1))
# lines(x=x$mids, y=(x$density*sum(x$counts)), col="red")
# lines(smooth.spline(x=x$mids, y=(x$density*sum(x$counts)), spar=0.35))
# dataTmp <- data.frame(x=x$mids,  y=(x$density*sum(x$counts)))
# lines(loess(y~x, data=dataTmp), col='red', lwd=2)

p <- gghistogram(aggSize_clust_dt, x = "clustSize", fill = "lightgray", binwidth=1, boundary=0,
            add = "mean", rug = F) +
  ggtitle("Distribution CTCF cluster size", subtitle=paste0("-d = ", d_clustSize)) +
  labs(x="cluster size", y="count")+
  scale_y_continuous(breaks=scales::pretty_breaks(n = 10)) + 
  # scale_x_continuous(breaks=scales::pretty_breaks(n = 10)) 
  scale_x_continuous(breaks=sort(unique(aggSize_clust_dt$clustSize)),
                     labels=sort(unique(aggSize_clust_dt$clustSize)))+
# scale_x_continuous(breaks=1:10,labels=sort(unique(aggSize_clust_dt$clustSize))) + 
  theme(
    plot.subtitle = element_text(size=14, face="italic", hjust=0.5),
    plot.title = element_text(size=16, face="bold", hjust=0.5),
    axis.text = element_text(size=12)
  )
#> Warning: Using `bins = 30` by default. Pick better value with the argument `bins`.

outFile <- file.path(outFolder, paste0("ctcf_clustSize_dist_d-", d_clustSize, ".svg"))
ggsave(p, filename = outFile, height=5, width=6)
cat(paste0("... written: ", outFile, "\n"))
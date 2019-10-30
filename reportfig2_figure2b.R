options(scipen=100)

# Rscript reportfig2_figure2b.R

script_name <- "reportfig2_figure2b.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(ggplot2)
require(reshape2)

registerDoMC(4)


plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <-  ifelse(plotType=="png", 600, 10)
axisCex <- 1.4

outFolder <- file.path("REPORTFIG2_FIGURE2b")
dir.create(outFolder, recursive = TRUE)

all_dt <- get(load(file.path("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")  ))
gene_signifThresh <- 0.05
tad_signifThresh <- 0.01

totNbr <- aggregate(entrezID ~hicds + exprds, data=all_dt, FUN=length)
colnames(totNbr)[colnames(totNbr) == "entrezID"] <- "totGenes"

tad_signifNbr <- aggregate(tad_adjCombPval ~hicds + exprds, data=all_dt, FUN=function(x) sum(x <= tad_signifThresh))
colnames(tad_signifNbr)[colnames(tad_signifNbr) == "tad_adjCombPval"] <- "nbr_tadSignif"

gene_signifNbr <- aggregate(adj.P.Val ~hicds + exprds, data=all_dt, FUN=function(x) sum(x <= gene_signifThresh))
colnames(gene_signifNbr)[colnames(gene_signifNbr) == "adj.P.Val"] <- "nbr_geneSignif"

plot_dt <- merge(gene_signifNbr, merge(totNbr, tad_signifNbr, by=c("hicds", "exprds")), by=c("hicds", "exprds"))

plot_dt_m <- melt(plot_dt, id=c("hicds", "exprds"))

plot_dt_m$value_log10 <- log10(plot_dt_m$value)

plot_dt_m$variable <- factor(as.character(plot_dt_m$variable), levels=rev(c("nbr_tadSignif", "nbr_geneSignif", "totGenes") ))
stopifnot(!is.na(plot_dt_m$variable))


outFile <- file.path(outFolder, paste0("nbrGenes_tot_signif_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))


boxplot(value_log10~variable, data=plot_dt_m,
        main=paste0("# of genes"),
        xlab="", 
        ylab="# of genes (log10)",
        names=c("total", "signif. gene level", "signif. TAD level" ),
        sub=paste0("adj. p-val gene level <= ", gene_signifThresh, "\n", "adj. p-val TAD level <= ", tad_signifThresh),
        cex.axis =axisCex,
        cex.lab =axisCex
)
mtext(side=3, text=paste0("all datasets - n=", length(unique(paste0(plot_dt$hicds, plot_dt$exprds)))))

legend("topright", legend= c(paste0("mean tot. = ", round(mean(plot_dt$totGenes),2), "\n", 
                                    "mean signif. gene level = ", round(mean(plot_dt$nbr_geneSignif),2),"\n",
                                    "mean signif. TAD level = ", round(mean(plot_dt$nbr_tadSignif),2))), bty="n")


foo <- dev.off()


cat(paste0("... written: ", outFile, "\n"))

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


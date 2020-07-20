# Rscript cmp_purity_data.R

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- file.path("CMP_PURITY_DATA")
dir.create(outFolder)

est_dt <- get(load("ALLTADS_AND_PURITY/log10/all_ds_corrPurity_dt.Rdata"))
colnames(est_dt)[colnames(est_dt)=="purityCorr"] <- "purityCorrEST"
agg_est_dt <- aggregate(purityCorrEST~dataset+region, data=est_dt, FUN=mean)

cpe_dt <- get(load("ALLTADS_AND_PURITY/CPE/log10/all_ds_corrPurity_dt.Rdata"))
colnames(cpe_dt)[colnames(cpe_dt)=="purityCorr"] <- "purityCorrCPE"
agg_cpe_dt <- aggregate(purityCorrCPE~dataset+region, data=cpe_dt, FUN=mean)

epic_dt <- get(load("ALLTADS_AND_PURITY/EPIC/log10/all_ds_corrPurity_dt.Rdata"))
colnames(epic_dt)[colnames(epic_dt)=="purityCorr"] <- "purityCorrEPIC"
agg_epic_dt <- aggregate(purityCorrEPIC~dataset+region, data=epic_dt, FUN=mean)

merge_cols <- c("dataset", "region", "entrezID")

merge_dt <- merge(epic_dt[,c(merge_cols, "purityCorrEPIC")], merge(est_dt[,c(merge_cols, "purityCorrEST")], cpe_dt[,c(merge_cols, "purityCorrCPE")], by=merge_cols), by=merge_cols)


merge_cols <- c("dataset", "region")
agg_merge_dt <- merge(agg_epic_dt[,c(merge_cols, "purityCorrEPIC")], 
                      merge(agg_est_dt[,c(merge_cols, "purityCorrEST")], agg_cpe_dt[,c(merge_cols, "purityCorrCPE")], by=merge_cols), by=merge_cols)


my_x <- merge_dt$purityCorrEST
my_y <- merge_dt$purityCorrCPE

outFile <- file.path(outFolder, "aranEST_aranCPE_geneLevel.png")
png(outFile, width=400, height=400)
plot(x=my_x,
     y=my_y, 
     pch=16,
     cex=0.7,
cex.main=1.2,
cex.lab=1.2,
cex.axis=1.2,
     xlab="Aran - ESTIMATES",
     ylab="Aran - CPE",
     main="Comparison purity data"
     )
mtext(side=3, text="gene level")
addCorr(x=my_y,y=my_y,legPos="topleft",bty="n")
curve(1*x, add=T,col="grey")
foo <- dev.off()

my_x <- merge_dt$purityCorrEST
my_y <- merge_dt$purityCorrEPIC


outFile <- file.path(outFolder, "aranEST_EPIC_geneLevel.png")
png(outFile, width=400, height=400)
plot(x=my_x,
     y=my_y, 
     pch=16,
     cex=0.7,
cex.main=1.2,
cex.lab=1.2,
cex.axis=1.2,
     xlab="Aran - ESTIMATES",
     ylab="EPIC",
     main="Comparison purity data"
)
mtext(side=3, text="gene level")
addCorr(x=my_x,y=my_y,legPos="topleft",bty="n")
curve(1*x, add=T,col="grey")
foo <- dev.off()

my_x <- merge_dt$purityCorrCPE
my_y <- merge_dt$purityCorrEPIC

outFile <- file.path(outFolder, "aranCPE_EPIC_geneLevel.png")
png(outFile, width=400, height=400)
plot(x=my_x,
     y=my_y,
     pch=16,
     cex=0.7,
cex.main=1.2,
cex.lab=1.2,
cex.axis=1.2,
     xlab="Aran - CPE",
     ylab="EPIC",
     main="Comparison purity data"
)
mtext(side=3, text="gene level")
addCorr(x=my_x,y=my_y,legPos="topleft",bty="n")
curve(1*x, add=T,col="grey")
foo <- dev.off()


###############

my_x <- agg_merge_dt$purityCorrEST
my_y <- agg_merge_dt$purityCorrCPE

outFile <- file.path(outFolder, "aranEST_aranCPE_tadLevel.png")
png(outFile, width=400, height=400)
plot(x=my_x,
     y=my_y,
     pch=16,
     cex=0.7,
cex.main=1.2,
cex.lab=1.2,
cex.axis=1.2,
     xlab="Aran - ESTIMATES",
     ylab="Aran - CPE",
     main="Comparison purity data"
)
mtext(side=3, text="TAD level")
addCorr(x=my_x,y=my_y,legPos="topleft",bty="n")
curve(1*x, add=T,col="grey")
foo <- dev.off()

my_x <- agg_merge_dt$purityCorrEST
my_y <- agg_merge_dt$purityCorrEPIC

outFile <- file.path(outFolder, "aranEST_EPIC_tadLevel.png")
png(outFile, width=400, height=400)
plot(x=my_x,
     y=my_y, 
     pch=16,
     cex=0.7,
cex.main=1.2,
cex.lab=1.2,
cex.axis=1.2,
     xlab="Aran - ESTIMATES",
     ylab="EPIC",
     main="Comparison purity data"
)
mtext(side=3, text="TAD level")
addCorr(x=my_x,y=my_y,legPos="topleft",bty="n")
curve(1*x, add=T,col="grey")
foo <- dev.off()

my_x <- agg_merge_dt$purityCorrCPE
my_y <- agg_merge_dt$purityCorrEPIC


outFile <- file.path(outFolder, "aranCPE_EPIC_tadLevel.png")
png(outFile, width=400, height=400)
plot(x=my_x,
     y=my_y,
     pch=16,
     cex=0.7,
cex.main=1.2,
cex.lab=1.2,
cex.axis=1.2,
     xlab="Aran - CPE",
     ylab="EPIC",
     main="Comparison purity data"
)
mtext(side=3, text="TAD level")
addCorr(x=my_x,y=my_y,legPos="topleft",bty="n")
curve(1*x, add=T,col="grey")
foo <- dev.off()








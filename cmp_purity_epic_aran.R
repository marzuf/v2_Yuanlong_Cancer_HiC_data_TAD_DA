aran_purity_file <- file.path("tcga_purity_aran2015.csv")
aran_purity_dt <- read.delim(aran_purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
pm <- purity_metrics[1]
aran_purity_dt$Sample.ID <- gsub("A$", "", aran_purity_dt$Sample.ID)

epic_purity_file <- file.path("EPIC_PURITY/all_epic_purity_data.Rdata")
epic_purity_data <- get(load(epic_purity_file))

epic_purity_dt <- as.data.frame(do.call(rbind, c(lapply(epic_purity_data, function(x) x[["infiltration_fraction"]]))))
                                                 
epic_purity_dt$Sample.ID <- rownames(epic_purity_dt)

stopifnot(any(epic_purity_dt$Sample.ID %in% aran_purity_dt$Sample.ID))

cmp_dt <- merge(aran_purity_dt, epic_purity_dt, by="Sample.ID", all.x=FALSE, all.y=FALSE)

nrow(aran_purity_dt)
nrow(epic_purity_dt)

nrow(cmp_dt)

outFolder <- "CMP_PURITY_EPIC_ARAN"
dir.create(outFolder, recursive = TRUE)

myHeight <- 400
myWidth <- 400
plotType <- "png"

outFile <- file.path(outFolder, paste0("epicOtherCells_vs_aranESTIMATE.", plotType))
do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
plot(
  otherCells ~ ESTIMATE,
  data = cmp_dt,
  pch=16,
  cex=0.7
)
legend(
  "topleft",
  legend=paste0("Pearson's coeff=", round(cor(cmp_dt$ESTIMATE, cmp_dt$otherCells, use="complete.obs"),4)),
  bty="n"
)

mtext(side=3, text = paste0("n=", nrow(cmp_dt)), cex=1.4)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))    

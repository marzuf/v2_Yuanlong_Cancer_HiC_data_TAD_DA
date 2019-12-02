aran_purity_file <- file.path("tcga_purity_aran2015.csv")
aran_purity_dt <- read.delim(aran_purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
pm <- purity_metrics[1]
aran_purity_dt$Sample.ID <- gsub("A$", "", aran_purity_dt$Sample.ID)

epic_purity_file <- file.path("EPIC_PURITY/all_epic_purity_data.Rdata")
epic_purity_data <- get(load(epic_purity_file))

samples_dt <- do.call(rbind, lapply(1:length(epic_purity_data), function(x) {
  tmp <- rownames(epic_purity_data[[x]][["infiltration_fraction"]])
  data.frame(
    hicds = dirname(names(epic_purity_data)[x]),
    exprds = basename(names(epic_purity_data)[x]),
    Sample.ID = tmp,
    stringsAsFactors = FALSE
  )
}))

epic_purity_dt <- as.data.frame(do.call(rbind, c(lapply(epic_purity_data, function(x) x[["infiltration_fraction"]]))))
                                                 
epic_purity_dt$Sample.ID <- rownames(epic_purity_dt)

stopifnot(any(epic_purity_dt$Sample.ID %in% aran_purity_dt$Sample.ID))

cmp_dt <- merge(aran_purity_dt, epic_purity_dt, by="Sample.ID", all.x=FALSE, all.y=FALSE)

nrow(aran_purity_dt)
nrow(epic_purity_dt)

nrow(cmp_dt)

source("subtype_cols.R")

stopifnot(cmp_dt$Sample.ID %in% samples_dt$Sample.ID)

cmp_samples_dt <- merge(cmp_dt, samples_dt, by=c("Sample.ID"))
cmp_samples_dt$cmp_type <- all_cmps[paste0(cmp_samples_dt$exprds)]
stopifnot(!is.na(cmp_samples_dt$cmp_type))
cmp_samples_dt$cmp_type_col <- ifelse(cmp_samples_dt$cmp_type == "norm_vs_tumor", "blue",
                                      ifelse(cmp_samples_dt$cmp_type == "subtypes", "green",
                                             ifelse(cmp_samples_dt$cmp_type == "wt_vs_mut", "red",NA)))
stopifnot(!is.na(cmp_samples_dt$cmp_type_col))

outFolder <- "CMP_PURITY_EPIC_ARAN"
dir.create(outFolder, recursive = TRUE)

myHeight <- 400
myWidth <- 400
plotType <- "png"

outFile <- file.path(outFolder, paste0("epicOtherCells_vs_aranESTIMATE.", plotType))
do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
plot(
  otherCells ~ ESTIMATE,
  data = cmp_samples_dt,
  pch=16,
  cex=0.7,
  cex.lab=1.4,
  cex.axis=1.4,
  col = cmp_samples_dt$cmp_type_col
)
legend(
  "topleft",
  legend=paste0("Pearson's coeff=", round(cor(cmp_dt$ESTIMATE, cmp_dt$otherCells, use="complete.obs"),4)),
  bty="n"
)
legend(
  "bottomleft",
  legend=c("norm_vs_tumor", "subtypes", "wt_vs_mut"),
  col = c("blue", "green", "red"),
  pch=16,
  bty="n"
)

mtext(side=3, text = paste0("n=", nrow(cmp_dt)), cex=1.4)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))    


cmp_samples_dt$corr_diff <- cmp_samples_dt$ESTIMATE - cmp_samples_dt$otherCells
cmp_samples_dt <- cmp_samples_dt[order(cmp_samples_dt$corr_diff, decreasing = TRUE),]
head(cmp_samples_dt[,c("hicds", "exprds", "otherCells", "ESTIMATE", "corr_diff")])

dis_dt <- cmp_samples_dt[cmp_samples_dt$ESTIMATE > 0.4 & cmp_samples_dt$otherCells < 0.4,]
dis_dt <- dis_dt[,c("hicds", "exprds", "otherCells", "ESTIMATE", "corr_diff")]
dis_dt <- na.omit(dis_dt)
dis_dt <- dis_dt[dis_dt$corr_diff >= 0.5,]
View(dis_dt[,c("hicds", "exprds", "otherCells", "ESTIMATE", "corr_diff")])
unique(dis_dt$exprds)
unique(dis_dt$exprds)
# [1] "TCGAlihc_norm_lihc"            "TCGAlihc_wt_mutCTNNB1"         "TCGAlgg_IDHwt_IDHmutnc"       
# [4] "TCGAgbm_classical_neural"      "TCGAgbm_classical_mesenchymal" "TCGAgbm_classical_proneural"  
# [7] "TCGAbrca_lum_bas"              "TCGAskcm_wt_mutBRAF"           "TCGAskcm_wt_mutCTNNB1"



dt1 <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
dt1 <- dt1[dt1$adjPvalComb <= 0.01,]
dt1 <- dt1[order(dt1$hicds, dt1$exprds, dt1$adjPvalComb),]
dt1$id <- file.path(dt1$hicds, dt1$exprds, dt1$region)
dt2 <- get(load("CREATE_FINAL_TABLE_PARTIAL/all_result_dt.Rdata"))
dt2 <- dt2[dt2$adjPvalComb <= 0.01,]
dt2 <- dt2[order(dt2$hicds, dt2$exprds, dt2$adjPvalComb),]
dt2$id <- file.path(dt2$hicds, dt2$exprds, dt2$region)

nrow(dt1)
nrow(dt2)
length(intersect(dt1$id, dt2$id))

dt1_10 <- do.call(rbind, by(dt1, list(dt1$hicds, dt1$exprds), FUN=head, n=10))
dt2_10 <- do.call(rbind, by(dt2, list(dt2$hicds, dt2$exprds), FUN=head, n=10))
nrow(dt1_10)
nrow(dt2_10)
length(intersect(dt1_10$id, dt2_10$id))

dt2 <- get(load("EPIC_PURITY/all_purity_values.Rdata"))
dt2 <- unlist(dt2, use.names = TRUE)
dt2_dt <- data.frame(
  id_code = gsub(".+TCGA-.+-.+-(.+)$", "\\1", names(dt2)),
  purity = as.numeric(dt2),
  stringsAsFactors = FALSE
)




boxplot(purity~id_code, data=dt2_dt)
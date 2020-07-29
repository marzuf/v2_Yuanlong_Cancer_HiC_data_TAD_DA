purity_file <- file.path("tcga_purity_aran2015.csv")
purity_dt <- read.delim(purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
pm <- purity_metrics[1]
purity_plot_name <- "aran"

est_dt <- get(load("PURITY_AVAILABLE_SAMPLES_PKG//ESTIMATE/all_ds_avPurity_dt.Rdata"))
cpe_dt <- get(load("PURITY_AVAILABLE_SAMPLES/CPE/all_ds_avPurity_dt.Rdata"))
stopifnot(setequal(est_dt$id_samp, cpe_dt$id_samp))

purity_dt <- purity_dt[,c("Sample.ID", pm)]
purity_dt <- na.omit(purity_dt)
nrow(purity_dt) 
# 7737 
purity_dt$Sample.ID_short <- substr(purity_dt$Sample.ID, start=1, stop=15)
cat(paste0(sum(unique(cpe_dt$id_samp) %in% purity_dt$Sample.ID_short), "/", length(unique(cpe_dt$id_samp)), "\n"))
#3519/3765

pkg_dt <- get(load("Tumor.purity.Rdata"))
pkg_dt <- pkg_dt[,c("Sample.ID", pm)]
pkg_dt <- na.omit(pkg_dt)
nrow(pkg_dt)
# 9364 
pkg_dt$Sample.ID_short <- substr(pkg_dt$Sample.ID, start=1, stop=15)
cat(paste0(sum(unique(cpe_dt$id_samp) %in% pkg_dt$Sample.ID_short), "/", length(unique(cpe_dt$id_samp)), "\n"))
#3523/3765

 # => these are exactly the same data !

pkg_est <- as.numeric(as.character(gsub(",", ".", pkg_dt$ESTIMATE)))
aran_est <- purity_dt$ESTIMATE
stopifnot(pkg_est == aran_est)

purity_file <- file.path("tcga_purity_aran2015.csv")
purity_dt <- read.delim(purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
purity_dt$Sample.vial <- substr(purity_dt$Sample.ID, start=16, stop=16)
purity_dt$Sample.short <- substr(purity_dt$Sample.ID, start=1, stop=15)
any(duplicated(purity_dt$Sample.short))

purity_dt <- purity_dt[,c("Sample.short", "Sample.vial", pm)]
purity_dt <- na.omit(purity_dt)

dt_a <- purity_dt[purity_dt$Sample.vial == "A",c("Sample.short", pm)]
dt_b <- purity_dt[purity_dt$Sample.vial == "B",c("Sample.short", pm)]
dt_c <- purity_dt[purity_dt$Sample.vial == "C",c("Sample.short", pm)]

dt_ab <- merge(dt_a, dt_b, by="Sample.short", all=FALSE, suffixes = c(".A", ".B")) # 30
dt_ac <- merge(dt_a, dt_c, by="Sample.short", all=FALSE, suffixes = c(".A", ".C")) # 4
dt_bc <- merge(dt_b, dt_c, by="Sample.short", all=FALSE, suffixes = c(".B", ".C")) # 0
plot(dt_ab$CPE.A, dt_ab$CPE.B, pch=16, cex=0.7)
plot(dt_ac$CPE.A, dt_ac$CPE.C, pch=16, cex=0.7)

any(dt_a$Sample.short %in% dt_b$Sample.short)
# [1] FALSE
any(dt_b$Sample.short %in% dt_a$Sample.short)
# [1] FALSE

any(dt_a$Sample.short %in% dt_c$Sample.short)
any(dt_b$Sample.short %in% dt_c$Sample.short)
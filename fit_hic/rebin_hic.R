#Rscript rebin_hic 10000 ENCSR444WCZ_A549_mat_chr10_10kb_ob.txt 20000 ENCSR444WCZ_A549_mat_chr10_20kb_ob.txt

infile="ENCSR444WCZ_A549_mat_chr10_10kb_ob.txt"

args <- commandArgs(trailingOnly = TRUE)
inBin <- as.numeric(args[1])
infile <- args[2]
outBin <- as.numeric(args[3])
outfile <- args[4]


dir.create(dirname(outfile), recursive = TRUE)

stopifnot(file.exists(infile))
stopifnot(is.numeric(inBin))
stopifnot(is.numeric(outBin))

in_dt <- read.delim(infile, header=F, col.names = c("binA", "binB", "count"), stringsAsFactors = FALSE)
stopifnot(in_dt$binA %% inBin == 0)
stopifnot(in_dt$binB %% inBin == 0)
in_dt$newA <- ceiling(in_dt$binA/outBin) * outBin
in_dt$newB <- ceiling(in_dt$binB/outBin) * outBin

# > head(in_dt)
# binA   binB     count   newA   newB
# 1  70000  70000 187.24530  80000  80000
# 2  80000  90000  40.83534  80000 100000
# 3  90000  90000 109.99067 100000 100000
# 4  80000 100000  25.72629  80000 100000
# 5  90000 100000  21.88238 100000 100000

out_dt <- aggregate(count~newA+newB, data=in_dt[,c("newA", "newB", "count")], FUN=sum)
stopifnot(is.numeric(out_dt$newA))
stopifnot(is.numeric(out_dt$newB))
stopifnot(out_dt$newB >= out_dt$newA)

cat(paste0("nrow(in_dt)\t", nrow(in_dt), "\n"))
cat(paste0("nrow(out_dt)\t", nrow(out_dt), "\n"))

write.table(out_dt, file=outfile, quote=F, col.names = F, row.names = F, sep="\t")
cat(paste0("... written: ", outfile, "\n"))



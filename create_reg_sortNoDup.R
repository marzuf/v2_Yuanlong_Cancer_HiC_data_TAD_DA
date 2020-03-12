startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
require(doMC)
registerDoMC(40)

# Rscript create_reg_sortNoDup.R 


outFolder <- file.path(paste0("CREATE_REG_SORTNODUP"))
dir.create(outFolder, recursive = TRUE)

reg_dt <- get(load("LOOK_VIPER/reg_dt.Rdata"))


stopifnot("regEntrezID" %in% colnames(reg_dt))
stopifnot("targetEntrezID" %in% colnames(reg_dt))


out_dt <- do.call(rbind, 
                 
                 by(reg_dt, reg_dt$regEntrezID, function(x) {
                   reg_entrez <- unique(x$targetEntrezID)
                   tf <- unique(x$regEntrezID)
                   stopifnot(length(tf) == 1)
                   all_cmbs <- combn(reg_entrez, m = 2)
                   stopifnot(nrow(all_cmbs) == 2)
                   stopifnot(ncol(all_cmbs) == 0.5*(length(reg_entrez)-1) * length(reg_entrez))
                   
                   data.frame(
                     gene1_tmp = all_cmbs[1,], 
                     gene2_tmp = all_cmbs[2,],
                    reg = tf,
                     stringsAsFactors = FALSE
                   )
                   
                 })
)

head(out_dt)

out_dt$gene1_tmp <- as.character(out_dt$gene1_tmp)
out_dt$gene2_tmp <- as.character(out_dt$gene2_tmp)
out_dt$gene1 <-as.character(pmin(out_dt$gene1_tmp, out_dt$gene2_tmp))
out_dt$gene2 <-as.character(pmax(out_dt$gene1_tmp, out_dt$gene2_tmp))

all_reg_pairs <- out_dt[,c("gene1", "gene2", "reg")]

stopifnot(all_reg_pairs$gene1 < all_reg_pairs$gene2)


outFile <- file.path(outFolder, paste0("all_reg_pairs.Rdata"))
save(all_reg_pairs, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))



######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
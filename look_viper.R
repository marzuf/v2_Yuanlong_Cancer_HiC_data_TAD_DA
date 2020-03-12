# Rscript look_viper.R

library(foreach)
library(doMC)

registerDoMC(40)

outFolder <- "LOOK_VIPER"
dir.create(outFolder)

# mrs_v0 <- get(load("TRIAL_VIPER/luad_mrs.Rdata"))
# mrs_v2 <- get(load("TRIAL_VIPER2/luad_mrs.Rdata"))
# 
# target_luad_mrs_v0 <- ledge(mrs_v0)
# target_luad_mrs_v2 <- ledge(mrs_v2)

# p.value FDR

library(aracne.networks)

viper_res <- read.delim("TRIAL_VIPER//mrs_res_table.txt", stringsAsFactors=FALSE, header=TRUE)
nrow(viper_res)
# 5057

pvalSignif_reg <- viper_res$Regulon[viper_res$p.value <= 0.05]
length(pvalSignif_reg)
# 1278
fdrSignif_reg <- viper_res$Regulon[viper_res$FDR <= 0.05]
length(fdrSignif_reg)
# 661

data(regulonluad)
length(regulonluad)

stopifnot(pvalSignif_reg %in% names(regulonluad))


reg_dt <- foreach(reg = pvalSignif_reg, .combine='rbind') %dopar% {
  
  all_targets <- regulonluad[[paste0(reg)]][["tfmode"]]
  stopifnot(length(all_targets) > 0)
  # take only the positive
  reg_targets <- names(all_targets)[all_targets > 0]
  if(length(reg_targets) == 0) return(NULL)
  
  data.frame(
    regEntrezID = as.character(reg),
    targetEntrezID = as.character(reg_targets),
    stringsAsFactors = FALSE
  )
  
}

outFile <- file.path(outFolder, "reg_dt.Rdata")
save(reg_dt, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


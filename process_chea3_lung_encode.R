chea3_lung_dt <- read.delim("chea3_lung_TFs_processed.txt", stringsAsFactors = FALSE, header=TRUE)
encode_lung_dt <- read.delim("chea3_ENCODE_processed.txt",  stringsAsFactors = FALSE, header=TRUE)

encode_lung_dt$regSymbol_2 <- gsub("(.+?)_.+", "\\1", encode_lung_dt$regSymbol)

all(encode_lung_dt$regSymbol %in% chea3_lung_dt$regSymbol)

sum(unique(chea3_lung_dt$regSymbol) %in% encode_lung_dt$regSymbol) / length(unique(chea3_lung_dt$regSymbol))

sum(unique(chea3_lung_dt$regSymbol) %in% encode_lung_dt$regSymbol_2) / length(unique(chea3_lung_dt$regSymbol))
sum(unique(chea3_lung_dt$regSymbol) %in% encode_lung_dt$regSymbol_2) 
length(unique(chea3_lung_dt$regSymbol))

# > sum(unique(chea3_lung_dt$regSymbol) %in% encode_lung_dt$regSymbol_2)
# [1] 118       ==> only 118 ??
# > length(unique(chea3_lung_dt$regSymbol))
# [1] 1620
permutFile <- "Panc1_rep12_40kb_TCGApaad_wt_mutKRAS_permutationsDT_1-100.Rdata"


hicds <- "Panc1_rep12_40kb"
exprds <- "TCGApaad_wt_mutKRAS"


script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"

pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))


stopifnot(file.exists(permutFile))
permutDT <- get(load(permutFile))

geneListFile <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
stopifnot(file.exists(geneListFile))
geneList <- get(load(geneListFile))

topTableFile <- file.path(pipFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")
stopifnot(file.exists(topTableFile))
topTable_DT <- get(load(topTableFile))

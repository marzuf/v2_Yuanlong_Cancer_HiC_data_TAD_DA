options(scipen=100)

setDir=""

# Rscript cat_tad_coord.R


hicds="ENCSR504OTV_transverse_colon_40kb"
exprds="TCGAcoad_msi_mss"

script_name <- "cat_tad_coord.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

outFolder <- "IGV_PLOTS"
dir.create(outFolder, recursive = TRUE)

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

pipFolder <- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

pvalThresh <- 0.01
matching_ratio_thresh <- 0.8

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2 | length(args) == 0)
hicds <- args[1]
exprds <- args[2]

all_hicds <- c("ENCSR504OTV_transverse_colon_40kb", "Rao_HCT-116_2017_40kb", "GSE105318_DLD1_40kb")
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))

final_table_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT[, paste0("signifEmpPval_", pvalThresh)] <- final_table_DT$adjPvalComb <= pvalThresh
final_table_DT <- final_table_DT[order(final_table_DT$adjPvalComb),]


snapshotDirectory <- "/home/marie/Desktop/IGV_plots_03.09/IGV_plots"



for(hicds in all_hicds)  {
  for(exprds in all_exprds[[paste0(hicds)]])  {

    
    file1 <- file.path(outFolder, paste0("IGV_plots_h3k4me1_", hicds, "_", exprds, ".batch"))
    file2 <- file.path(outFolder, paste0("IGV_plots_h3k27ac_", hicds, "_", exprds, ".batch"))
    file.remove(file1)
    file.remove(file2)
    
    myconH3K4me1 <- file(file1, "a")
    myconH3K27ac <- file(file2, "a")
    
    snapshotDirectory <- "/home/marie/Desktop/IGV_plots_03.09/IGV_plots"
    
    writeLines(paste("snapshotDirectory", snapshotDirectory), con=myconH3K4me1)
    writeLines(paste("snapshotDirectory", snapshotDirectory), con=myconH3K27ac)
    writeLines(paste("maxPanelHeight 1200"), con=myconH3K4me1)
    writeLines(paste("maxPanelHeight 1200"), con=myconH3K27ac)
    # goto chr10:85879000-86321001
    # snapshot H3K4me1_chr10_85880001_86320000_slop1000.png
    
    
    tadFile <- file.path(pipFolder, hicds, "genes2tad", "all_assigned_regions.txt")
    stopifnot(file.exists(tadFile))
    tad_DT <- read.delim(tadFile, header=F, col.names = c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
    
    # chr10   85880000        86320001        chr10_TAD113
    
    
    saveDT <- final_table_DT[final_table_DT[, paste0("signifEmpPval_", pvalThresh)] & final_table_DT$hicds == hicds & final_table_DT$exprds == exprds, c("hicds", "exprds", "region", "meanLogFC", paste0("signifEmpPval_", pvalThresh)) ]
    outfile <- file.path(outFolder, paste0(hicds, "_", exprds, "_TAD_meanLogFC.txt"))
    write.table(saveDT, col.names = FALSE, row.names = TRUE, sep="\t", quote=F, append=F, file=outfile)
    
    signifTADs <- final_table_DT$region[final_table_DT[, paste0("signifEmpPval_", pvalThresh)] & final_table_DT$hicds == hicds & final_table_DT$exprds == exprds ]
    
    selectDT <- tad_DT[tad_DT$region %in% signifTADs, c("chromo", "start", "end", "region")]
    selectDT$start <- selectDT$start - 1
    selectDT$end <- selectDT$end + 1
    
    
    write.table(selectDT[order(factor(selectDT$region, levels=signifTADs)),], file=file.path(outFolder, paste0(hicds, "_", exprds, "_TADcoord.bed")), col.names = FALSE, row.names = FALSE, sep="\t", quote=F, append=F)

    for(histMark in c("H3K4me1", "H3K27ac")) {
      for(i in 1:nrow(selectDT)) {
        tad_id <- selectDT$region[i]
        tad_chromo <- selectDT$chromo[i]
        tad_start <- selectDT$start[i]
        tad_end <- selectDT$end[i]
        
        writeLines(paste0("goto ", tad_chromo, ":", tad_start, "-", tad_end), con=get(paste0("mycon", histMark)))
        writeLines(paste0("snapshot ", histMark, "_", hicds, "_", exprds, "_", tad_id, "_", tad_chromo, "_", tad_start, "_", tad_end, "_slop1000.png"), con=get(paste0("mycon", histMark)))
        
        
      }  #end-for tad
    }# end for histmark
  }#end-for exprds
}#end-for hicds

close(myconH3K4me1)
close(myconH3K27ac)

########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript dataset_purity_cmp.R\n"))

script_name <- "dataset_purity_cmp.R"

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript dataset_purity_cmp.R
# Rscript dataset_purity_cmp.R EPIC

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- TRUE


args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  purity_ds <- "" 
  file_suffix <- ""
} else {
  stopifnot(length(args) == 1)
  purity_ds <- args[1]
  file_suffix <- paste0("_", purity_ds)
}


outFolder <- file.path(paste0("DATASET_PURITY_CMP", file_suffix))
dir.create(outFolder, recursive = TRUE)

myHeightGG <- 7
myWidthGG <- 9
plotType <- "png"

mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")

if(purity_ds == "") {
  purity_file <- file.path("tcga_purity_aran2015.csv")
  purity_dt <- read.delim(purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
  pm <- purity_metrics[1]
  # all the ranks are between 1 and 0
  
} else if(purity_ds == "EPIC") {
  
  
  purity_file <- file.path("EPIC_PURITY/all_epic_purity_data.Rdata")
  epic_purity_data <- get(load(purity_file))
  purity_dt <- as.data.frame(do.call(rbind, c(lapply(epic_purity_data, function(x) x[["infiltration_fraction"]]))))
  purity_dt <- data.frame(purity_dt)
  purity_metrics <- colnames(purity_dt) #"Bcells"      "CAFs"        "CD4_Tcells"  "CD8_Tcells"  "Endothelial" "Macrophages" "NKcells"     "otherCells" 
  pm <- "otherCells"
  purity_dt$Sample.ID <- rownames(purity_dt)
  purity_dt$Sample.ID <- gsub("\\.", "-", purity_dt$Sample.ID)
  
  
  
} else{
  stop("---invalid DS\n")
}

all_hicds <- list.files(pipFolder)
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]

# all_ds=all_ds[1]
if(buildTable) {
  dataset_samp_cond_purity_dt <- foreach(ds = all_ds, .combine='rbind') %do% {
    hicds <- dirname(ds)
    exprds <- basename(ds)
    settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
    stopifnot(file.exists(settingFile))
    source(settingFile)
    samp1 <- get(load(file.path(setDir, sample1_file)))
    samp2 <- get(load(file.path(setDir, sample2_file)))
    pur_samp1 <- samp1[samp1 %in% purity_dt$Sample.ID | paste0(samp1, "A") %in% purity_dt$Sample.ID]
    cat(paste0("For ", cond1, " - available samples:\t", length(pur_samp1), "/", length(samp1), "\n"))
    pur_samp2 <- samp2[samp2 %in% purity_dt$Sample.ID | paste0(samp2, "A") %in% purity_dt$Sample.ID]
    cat(paste0("For ", cond2, " - available samples:\t", length(pur_samp2), "/", length(samp2), "\n"))
    
    if(length(pur_samp1) == 0 & length(pur_samp2) == 0) {
      return(NULL)
    }
    
    pur2_samp1 <- purity_dt$Sample.ID[purity_dt$Sample.ID %in% samp1  | purity_dt$Sample.ID %in% paste0(samp1, "A") ]
    stopifnot(length(pur2_samp1) == length(pur_samp1))
    pur2_samp2 <- purity_dt$Sample.ID[purity_dt$Sample.ID %in% samp2  | purity_dt$Sample.ID %in% paste0(samp2, "A") ]
    stopifnot(length(pur2_samp2) == length(pur_samp2))
    stopifnot(setequal(gsub("A$", "", pur2_samp1), pur_samp1))
    stopifnot(setequal(gsub("A$", "", pur2_samp2), pur_samp2))
    
    ds_pur_dt <- purity_dt[purity_dt$Sample.ID %in% c(pur2_samp1, pur2_samp2),]
    ds_pur_dt <- ds_pur_dt[,c("Sample.ID", purity_metrics)]
    ds_pur_dt_m <- melt(ds_pur_dt, id="Sample.ID")
    ds_pur_dt_m$cond <- ifelse(ds_pur_dt_m$Sample.ID %in% pur2_samp1, cond1,
                               ifelse(ds_pur_dt_m$Sample.ID %in% pur2_samp2, cond2, NA))
    stopifnot(!is.na(ds_pur_dt_m$cond))
    ds_pur_dt_m$hicds <- hicds
    ds_pur_dt_m$exprds <- exprds
    ds_pur_dt_m
  } # end-foreach iterating over datasets
  
  outFile <- file.path(outFolder, "dataset_samp_cond_purity_dt.Rdata")
  save(dataset_samp_cond_purity_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
} else {
  outFile <- file.path(outFolder, "dataset_samp_cond_purity_dt.Rdata")
  dataset_samp_cond_purity_dt <- get(load(outFile))
  
}

pm = purity_metrics[1]
for(pm in purity_metrics) {
  
  plot_dt <- dataset_samp_cond_purity_dt[dataset_samp_cond_purity_dt$variable == pm,]
  plot_dt$dataset <- paste0(plot_dt$hicds, "\n", plot_dt$exprds)
  plot_dt$subTit <- paste0(plot_dt$hicds, "-", plot_dt$exprds)
  ds = unique(plot_dt$dataset)[1]
  for(ds in unique(plot_dt$dataset)) {
    
    plot_dt_ds <- plot_dt[plot_dt$dataset == ds,]
    
    subtit <- unique(plot_dt_ds$subTit)
    
    p_var <-  ggplot(plot_dt_ds, aes(x = cond, y = value, fill = cond)) +
      geom_boxplot()+
      facet_grid(~dataset, switch="x") +
      coord_cartesian(expand = FALSE) +
      ggtitle(paste0(pm), subtitle = paste0(subtit))+
      scale_x_discrete(name="")+
      scale_y_continuous(name=paste0(),
                         breaks = scales::pretty_breaks(n = 20))+
      scale_fill_brewer(palette="Set1")+
      labs(fill  = "") +
      theme( # Increase size of axis lines
        strip.text = element_text(size = 12),
        # top, right, bottom and left
        # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", size=16),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.minor.y = element_line(colour = "grey"),
        strip.text.x = element_text(size = 10),
        axis.line.x = element_line(size = .2, color = "black"),
        axis.line.y = element_line(size = .3, color = "black"),
        axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        legend.background =  element_rect(),
        legend.key = element_blank(),
        legend.title = element_text(face="bold")
      )
    
    outFile <- file.path(outFolder, paste0(pm, "_", subtit, "_condComp_boxplot.", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
  } # end for iterating over the datasets
  
  # p_var <-  ggplot(plot_dt, aes(x = cond, y = value, fill = cond)) +
  #   geom_boxplot()+
  #   facet_grid(~dataset, switch="x") +
  #   coord_cartesian(expand = FALSE) +
  #   ggtitle(paste0(pm), subtitle = paste0(""))+
  #   scale_x_discrete(name="")+
  #   scale_y_continuous(name=paste0(),
  #                      breaks = scales::pretty_breaks(n = 20))+
  #   # scale_fill_brewer(palette="YlOrRd")+
  #   labs(fill  = "meanCorr type") +
  #   theme( # Increase size of axis lines
  #     strip.text = element_text(size = 12),
  #     # top, right, bottom and left
  #     # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
  #     plot.title = element_text(hjust = 0.5, face = "bold", size=16),
  #     plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
  #     panel.grid = element_blank(),
  #     panel.grid.major.y = element_line(colour = "grey"),
  #     panel.grid.minor.y = element_line(colour = "grey"),
  #     strip.text.x = element_text(size = 10),
  #     axis.line.x = element_line(size = .2, color = "black"),
  #     axis.line.y = element_line(size = .3, color = "black"),
  #     axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
  #     axis.text.x = element_blank(),
  #     axis.ticks.x = element_blank(),
  #     axis.title.y = element_text(color="black", size=12),
  #     axis.title.x = element_text(color="black", size=12),
  #     panel.border = element_blank(),
  #     panel.background = element_rect(fill = "transparent"),
  #     legend.background =  element_rect(),
  #     legend.key = element_blank(),
  #     legend.title = element_text(face="bold")
  #   )
  # 
  # outFile <- file.path(outFolder, paste0(pm, "_allDS_condComp_boxplot.", plotType))
  # ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG*3)
  # cat(paste0("... written: ", outFile, "\n"))
  
}  # end for iterating over purity metrics




##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))

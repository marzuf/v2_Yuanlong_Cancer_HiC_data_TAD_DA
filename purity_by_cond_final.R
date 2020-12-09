########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript allTADs_and_purity_final.R\n"))

script_name <- "purity_by_cond_final.R"

# _final -> discussion with Giovanni 04.08.2020 -> take Aran CPE data
# corrected compared to some of the previous versions -> if only non-"A" vial, take the "A" vial
# if multiple vials -> take "A" vials


purity_file <- file.path("tcga_purity_aran2015.csv")
purity_dt <- read.delim(purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
dt0 <- na.omit(purity_dt[,c("Sample.ID", "CPE")])
require(TCGAbiolinks)
dt1 <- na.omit(Tumor.purity[,c("Sample.ID", "CPE")])
dt1$CPE <- as.numeric(as.character(gsub(",", ".", as.character(dt1$CPE))))
dt1 <- na.omit(dt1)
dt0 <- dt0[order(dt0$Sample.ID),]
dt1 <- dt1[order(dt1$Sample.ID),]
dt1$Sample.ID <- as.character(dt1$Sample.ID)
stopifnot(all.equal(dt0, dt1))
# TRUE

curr_hicds <- "ENCSR489OCU_NCI-H460_40kb"
curr_exprds <- "TCGAluad_mutKRAS_mutEGFR"  


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_cols[all_cols=="green"] <- "darkgreen"

do_plot <- function(my_x, my_y, ...) {
  plot(x=my_x,
       y=my_y,
       pch=16,
       cex=0.7,
       cex.axis=1.2,
       cex.main=1.2,
       cex.lab=1.2,
       ...)
  addCorr(my_x, my_y, bty="n")
}

# Rscript purity_by_cond_final.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- F

fontFamily <- "Arial"

myHeight <- 400 
myWidth <- 400
plotType <- "png"
plotCex <- 1.4
mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")

myHeightGG <- 6
myWidthGG <- 9

# plotType <- "svg"
# myHeightGG <- 7
# myWidthGG <- 11
# plotCex <- 1.2

corMet <- "pearson"

transfExpr <- "log10"
logOff <- 0.01


# for plotting
signifThresh <- 0.01
highCorrThresh <- 0.6

# to quickly retrieve tad-level stat.
all_result_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  purity_ds <- args[1]  
} else{
  purity_ds <- ""
}

script0_name <- "0_prepGeneData"

purity_ds <- "aran"
file_suffix <- ""
purity_file <- file.path("tcga_purity_aran2015.csv")
purity_dt <- read.delim(purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
pm <- "CPE"
purity_plot_name <- "Aran - CPE"

# discussion with Giovanni 04.08.2020 -> if multiple vials, take the 1st one (e.g. if A and B, take A)
# -> remove duplicated substr sample.id will do the stuff
purity_dt <- purity_dt[,c("Sample.ID", pm)]
purity_dt <- na.omit(purity_dt)
purity_dt$Sample.ID <- as.character(purity_dt$Sample.ID)
purity_dt <- purity_dt[order(purity_dt$Sample.ID),]
purity_dt$Sample.ID_withVial <- purity_dt$Sample.ID
purity_dt$Sample.ID <- substr(start=1, stop=15, purity_dt$Sample.ID)
purity_dt_init <- purity_dt
purity_dt <- purity_dt[!duplicated(purity_dt$Sample.ID),]
discd_samps <- setdiff(purity_dt_init$Sample.ID, purity_dt$Sample.ID)
stopifnot(length(discd_samps) == 0 )
discd_vials <- setdiff(purity_dt_init$Sample.ID_withVial, purity_dt$Sample.ID_withVial)
cat(paste0("discd vials:\t", length(discd_vials), "\n"))
stopifnot(length(discd_vials) == nrow(purity_dt_init) - nrow(purity_dt))
stopifnot(!grepl("A$", discd_vials))
stopifnot(any(grepl("A$", purity_dt$Sample.ID_withVial)))
stopifnot(!duplicated(purity_dt$Sample.ID))
purity_dt$Sample.ID_withVial <- NULL

outFolder <- file.path("PURITY_BY_COND_FINAL", purity_ds, pm, transfExpr)
dir.create(outFolder, recursive = TRUE)

cat(paste0("!!! > purity metric: ", pm, "\n"))

all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]

cat(paste0(">>> purity metric\t=\t", pm, "\n"))

mainFolder <- "."

# all_ds=all_ds[1]

ex_hicds <- "ENCSR489OCU_NCI-H460_40kb"
ex_exprds <- "TCGAluad_norm_luad"
ex_TAD <- "chr11_TAD390"

# all_ds = file.path(ex_hicds, ex_exprds)

if(buildTable){
  
  all_ds_corrPurity_data <- foreach(ds = all_ds) %dopar% {
    
    cat(paste0("... start: ", ds, "\n"))
    
    hicds <- dirname(ds)
    exprds <- basename(ds)
    
    settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
    stopifnot(file.exists(settingFile))
    source(settingFile)
    
    samp1 <- get(load(file.path(setDir, sample1_file)))
    samp2 <- get(load(file.path(setDir, sample2_file)))
    
    pur_samp1 <- samp1[samp1 %in% purity_dt$Sample.ID ]
    cat(paste0("For ", cond1, " - available samples:\t", length(pur_samp1), "/", length(samp1), "\n"))
    
    pur_samp2 <- samp2[samp2 %in% purity_dt$Sample.ID ]
    cat(paste0("For ", cond2, " - available samples:\t", length(pur_samp2), "/", length(samp2), "\n"))
    
    if(length(pur_samp1) == 0 & length(pur_samp2) == 0) return(NULL)
    
    pur2_samp1 <- purity_dt$Sample.ID[purity_dt$Sample.ID %in% samp1 ]
    stopifnot(length(pur2_samp1) == length(pur_samp1))
    
    pur2_samp2 <- purity_dt$Sample.ID[purity_dt$Sample.ID %in% samp2 ]
    stopifnot(length(pur2_samp2) == length(pur_samp2))
    
    stopifnot(setequal(pur2_samp1, pur_samp1))
    stopifnot(setequal(pur2_samp2, pur_samp2))
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    
    geneList <- get(load(file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")))
	pipeline_geneList <- geneList
    stopifnot(geneList %in% g2t_dt$entrezID)
    g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
    g2t_dt$region <- as.character(g2t_dt$region)

    
    fpkm_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_fpkmDT.Rdata")
    stopifnot(file.exists(fpkm_file))
    fpkm_dt <- get(load(fpkm_file))  
    # DO THE SAME AS IN THE MANUSCRIPT_FIGURES SCRIPTS (cf discussion marco)
    # changed here 11.03.2020 -> should be normalized sample-wise
    fpkm_dt2 <- apply(fpkm_dt, 2, function(x)x/sum(x))
    # stopifnot(colSums(fpkm_dt2) == 1)
    stopifnot(abs(colSums(fpkm_dt2) - 1) <= 10^-4)
    # and then multiply by 10^6 to have FPKM
    fpkm_dt2 <- fpkm_dt2*10^6
    fpkm_dt2 <- data.frame(fpkm_dt2, check.names = FALSE)
    stopifnot(dim(fpkm_dt) == dim(fpkm_dt2))
    stopifnot(rownames(fpkm_dt) == rownames(fpkm_dt2))
    stopifnot(colnames(fpkm_dt) == colnames(fpkm_dt2))
    fpkm_dt <- fpkm_dt2
    
    stopifnot(names(geneList) %in% rownames(fpkm_dt))
    
    ### -> do the same as for the GIMAPs but for all TADs !
    
    purity_values <- setNames(purity_dt[purity_dt$Sample.ID %in% pur2_samp1 | purity_dt$Sample.ID %in% pur2_samp2,paste0(pm)],
                              purity_dt[purity_dt$Sample.ID %in% pur2_samp1 | purity_dt$Sample.ID %in% pur2_samp2,paste0("Sample.ID")])
    
    stopifnot( gsub("A$", "", names(purity_values)) == names(purity_values) )
    
    ds_purity_dt <- data.frame(
      dataset=ds,
      sampID = names(purity_values),
      purity = purity_values,
      stringsAsFactors = FALSE
    )
    rownames(ds_purity_dt) <- NULL
    
    ds_purity_dt$cond <- ifelse(ds_purity_dt$sampID %in% cond1_ID, cond1,
                             ifelse(ds_purity_dt$sampID %in% cond2_ID, cond2, NA))
    stopifnot(!is.na(ds_purity_dt$cond))
    
    tad = unique(g2t_dt$region)[1]
    
    all_tads_dt <- foreach(tad = unique(g2t_dt$region), .combine='rbind') %dopar% {
      
      tad_entrez <- g2t_dt$entrezID[g2t_dt$region == tad]
      stopifnot(tad_entrez %in% geneList)

      tad_entrez_fpkm <- names(geneList)[geneList %in% tad_entrez]
      stopifnot(tad_entrez_fpkm %in% rownames(fpkm_dt))
      tad_fpkm_dt <- fpkm_dt[tad_entrez_fpkm,]
      stopifnot(nrow(tad_fpkm_dt) == length(tad_entrez_fpkm))
      
      stopifnot(pur_samp1 %in% colnames(fpkm_dt)); stopifnot(pur_samp2 %in% colnames(fpkm_dt))
      stopifnot(pur2_samp1 %in% purity_dt$Sample.ID); stopifnot(pur2_samp2 %in% purity_dt$Sample.ID)
      
      all_samps <- sort(c(pur_samp1, pur_samp2))
      stopifnot(!duplicated(all_samps))
      all_samps2 <- sort(c(pur2_samp1, pur2_samp2))
      stopifnot(!duplicated(all_samps2))
      stopifnot(all_samps2 == all_samps)
      
      
      stopifnot(setequal(names(purity_values), all_samps2))
      stopifnot(setequal(names(purity_values), all_samps))
      
      
      stopifnot(all_samps %in% colnames(tad_fpkm_dt))
      
      tad_dt <- data.frame(t(tad_fpkm_dt[,all_samps]), check.names=FALSE)
      
      stopifnot(is.numeric(unlist(tad_dt)))
      
      if(!is.null(transfExpr)) {
        if(grepl("log", transfExpr)) {
          tad_dt_2 <- do.call(transfExpr, list(tad_dt+logOff))
          stopifnot(dim(tad_dt_2)==dim(tad_dt))
          stopifnot(rownames(tad_dt_2) == rownames(tad_dt))
          stopifnot(colnames(tad_dt_2) == colnames(tad_dt))
          tad_dt <- tad_dt_2
        } else {stop("unnknown\n")}
        labTransf <- transfExpr
      } else {
        labTransf <- ""
      }
      
      stopifnot(colnames(tad_dt) %in% names(pipeline_geneList))
      stopifnot(!duplicated(pipeline_geneList))
      # reback to gff dt entrezID -> needed to add column if missing gmaps
      colnames(tad_dt) <- pipeline_geneList[colnames(tad_dt)]
      stopifnot(colnames(tad_dt) %in% pipeline_geneList)
      
      stopifnot(setequal(colnames(tad_dt),  tad_entrez))
      tad_dt$sampID <- rownames(tad_dt)
      rownames(tad_dt) <- NULL
      
      purity_expr_dt <- merge(ds_purity_dt, tad_dt, by="sampID", all=TRUE)
      stopifnot(ncol(purity_expr_dt) == length(tad_entrez) + 4) ### added cond !
      purity_expr_dt <- purity_expr_dt[,c("dataset", "sampID", "purity", "cond", tad_entrez)]
      
      if(hicds==ex_hicds & exprds==ex_exprds & tad == ex_TAD) {
        for(i_g in tad_entrez) {
          outFile <- file.path(outFolder, paste0(ex_hicds, "_", ex_exprds, "_", ex_TAD, "_", i_g, "_expr_corr.", plotType))
          do.call(plotType, list(outFile, height=myHeight, width=myWidth))
          do_plot(my_x=purity_expr_dt[,"purity"],
                  my_y=purity_expr_dt[,paste0(i_g)],
                  main=paste0(hicds, " - ", exprds),
                  xlab="purity", ylab=paste0(i_g, " expr. ", labTransf))
          mtext(side=3, text=paste0(tad))
          foo <- dev.off()
          cat(paste0("... written: ", outFile, "\n"))
        }
      }
      # stop("ok \n")
      # correlation each column with purity 
      
      purity_expr_dt <- na.omit(purity_expr_dt)
      
      stopifnot(!duplicated(purity_expr_dt$sampID))
      
      if(nrow(purity_expr_dt) == 0) return(NULL)
      
      
      # all_purityCors <- apply(purity_expr_dt[,tad_entrez],2, function(col) cor(col, purity_expr_dt$purity, method=corMet))
      # 
      # stopifnot(!is.na(all_purityCors))
      # 
      # data.frame(
      #   dataset=ds,
      #   nSampWithPurity=length(purity_expr_dt$sampID),
      #   region = tad,
      #   entrezID= names(all_purityCors),
      #   purityCorr = as.numeric(all_purityCors),
      #   stringsAsFactors = FALSE
      # )
      
      ### CHANGED HERE 23.09.20 to have by condition
      
      out_dt <- do.call(rbind, by(purity_expr_dt,purity_expr_dt$cond, function(sub_dt){
        
        all_purityCors <- apply(sub_dt[,tad_entrez],2, function(col) cor(col, sub_dt$purity, method=corMet))
        
        stopifnot(!is.na(all_purityCors))
        data.frame(
          dataset=ds,
          cond = unique(sub_dt$cond),
          nSampWithPurity=length(sub_dt$sampID),
          region = tad,
          entrezID= names(all_purityCors),
          purityCorr = as.numeric(all_purityCors),
          stringsAsFactors = FALSE
        )
      } ))
      
      rownames(out_dt) <- NULL
      out_dt
    }
    list(
      ds_purity_dt=ds_purity_dt,
      all_tads_dt=all_tads_dt
    )
  }
  outFile <- file.path(outFolder,"all_ds_corrPurity_data.Rdata" )
  save(all_ds_corrPurity_data, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder,"all_ds_corrPurity_data.Rdata" )
  all_ds_corrPurity_data <- get(load(outFile))
}

all_dt <- do.call(rbind, lapply(all_ds_corrPurity_data, function(x)x[["ds_purity_dt"]]))
curr_dt <- all_dt[all_dt$dataset == file.path(curr_hicds, curr_exprds),]
curr_dt$cond_id <- curr_dt$cond
curr_dt$cond1 <- gsub("TCGA.+_(.+)_.+", "\\1", basename(curr_dt$dataset))
curr_dt$cond2 <- gsub("TCGA.+_.+_(.+)", "\\1", basename(curr_dt$dataset))
curr_dt$cond <- ifelse(curr_dt$cond_id == curr_dt$cond1, "cond1", 
                               ifelse(curr_dt$cond_id == curr_dt$cond2, "cond2", NA))
stopifnot(!is.na(curr_dt$cond))
curr_dt$cond_lab <- paste0(curr_dt$cond, " (", curr_dt$cond_id, ")")

plotTit <- paste0(curr_hicds, " - ", curr_exprds, ": sample purity level")
mySub <- paste0(paste0(names(table(curr_dt$cond_id)), ": ", as.numeric(table(curr_dt$cond_id))), collapse="; ")

require(ggpubr)

legTitle <- ""

save(curr_dt, file=file.path(outFolder, "curr_dt.Rdata"), version=2)
# stop("-ok\n")

p2 <- ggdensity(curr_dt,
                x = "purity",
                y = "..density..",
                # combine = TRUE,                  # Combine the 3 plots
                xlab = "sample purity", 
                # add = "median",                  # Add median line. 
                rug = FALSE,                      # Add marginal rug
                color = "cond_lab", 
                fill = "cond_lab",
                palette = "jco"
) +
  # scale_color_manual(values=my_cols)+
  # scale_fill_manual(values=my_cols)  +
  ggtitle(plotTit, subtitle = mySub)+
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") + 
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold")
  ) 


outFile <- file.path(outFolder, paste0(curr_hicds, "_", curr_exprds, "_cond1_cond2_density.", "svg"))
ggsave(p2, filename = outFile, height=5, width=7)


p2 <- ggboxplot(curr_dt,
                x = "cond_lab",
                y = "purity",
                # combine = TRUE,                  # Combine the 3 plots
                xlab = "sample purity", 
                # add = "median",                  # Add median line. 
                # rug = FALSE,                      # Add marginal rug
                color = "cond_lab", 
                fill = "cond_lab",
                palette = "jco"
) +
  # scale_color_manual(values=my_cols)+
  # scale_fill_manual(values=my_cols)  +
  ggtitle(plotTit, subtitle = mySub)+
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") + 
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold")
  ) 

outFile <- file.path(outFolder, paste0(curr_hicds, "_", curr_exprds, "_cond1_cond2_boxplot.", "svg"))
ggsave(p2, filename = outFile, height=5, width=7)



######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

# purity data form  Genomic Landscape of Non-Small Cell Lung Cancer in Smokers and Never-Smokers https://doi.org/10.1016/j.cell.2012.08.024
dataDT <-	read.table(textConnection("
Case	Gender	Smoking_Status	Dominant_Secondary_Clone_VAFs	Tumor_Purity__Dominant_Clone_VAFs	Tumor_Purity__X_VAFs	Clonality_Status
LUC1	male	light_smoker	12.7%	25.4%	30.7%	monoclonal
LUC2	female	smoker	22.5%/12.9%	45.0%	n/a	biclonal
LUC4	male	smoker	14.9%	29.8%	29.8%	monoclonal
LUC6	female	never-smoker	24.7%	49.4%	n/a	monoclonal
LUC7	female	never-smoker	21.3%10.8%	42.6%	n/a	biclonal
LUC8	female	smoker	22.7%	45.4%	n/a	monoclonal
LUC9	female	smoker	41.1%/20.4%	82.2%	n/a	biclonal
LUC10	male	smoker	43.1%/21.9%	86.2%	71.4%	biclonal
LUC11	male	never-smoker	28.8%	57.6%	41.0%	monoclonal
LUC12	male	smoker	10.4%	20.8%	24.3%	monoclonal
LUC13	male	smoker	41.9%/21.3%	83.8%	58.6%	biclonal
LUC14	female	smoker	29.5%/16.2%	59.0%	n/a	biclonal
LUC15	female	never-smoker	19.2%/10.8%	38.4%	n/a	biclonal
LUC16	female	never-smoker	47.2%/15.3%	94.4%	n/a	biclonal
LUC17	female	smoker	13.9%/11.2%	27.8%	n/a	biclonal
LUC18	male	smoker	18.8%/9.8%	37.6%	31.5%	biclonal
LUC20	female	smoker	39.9%	79.8%	n/a	monoclonal
                                    "),
                     header=TRUE)

dataDT$Tumor_Purity__Dominant_Clone_VAFs <- as.numeric(as.character(gsub("%", "", dataDT$Tumor_Purity__Dominant_Clone_VAFs)))
dataDT$Tumor_Purity__X_VAFs <- as.numeric(as.character(gsub("%", "", dataDT$Tumor_Purity__X_VAFs)))

aggregate(Tumor_Purity__Dominant_Clone_VAFs~Smoking_Status, data=dataDT, FUN=mean, na.rm=T)
aggregate(Tumor_Purity__X_VAFs~Smoking_Status, data=dataDT, FUN=mean, na.rm=T)



# Rscript meth_data_for_tad.R 

require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)

outFolder <- file.path("METH_DATA_FOR_TAD")
dir.create(outFolder, recursive = TRUE)

buildTable <- TRUE

all_hicds <- c("ENCSR444WCZ_A549_40kb", "LG1_40kb", "LG2_40kb", "ENCSR489OCU_NCI-H460_40kb")
hicds="LG1_40kb"
gene_oi <- "AKR1C1"

all_files <- list.files("encode_data", pattern="\\.bed$", full.names=TRUE)

runFolder <- ".."

normQt <- 0.95
qvalThresh <- 0.01
pvalThresh <- 0.01

plotType <-  "svg"
myHeight <- 7
myWidth <- 10
plotCex <- 1.4

# my_cols <- setNames(pal_jama()(5)[c(3, 2)], c(check_exprds, "other"))
# my_cols <- setNames(pal_jama()(5)[c(3, 2)], c(check_exprds, "other"))

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

stopifnot(gene_oi %in% entrez2symb)

entrez_oi <- names(entrez2symb)[entrez2symb == gene_oi]

encode_annot_dt <- read.delim("metadata_filter.csv", header=FALSE, sep=",", stringsAsFactors = FALSE)
bed2annot <- setNames(as.character(encode_annot_dt$V22), as.character(encode_annot_dt$V1))
bed2annot <- gsub("-human", "",  bed2annot)
bed2cl <- setNames(as.character(encode_annot_dt$V10), as.character(encode_annot_dt$V1))

cl_levels <- c("lung",
               "upper lobe of left lung",
               "AG04450", # fetal
               "IMR-90" , # normal
               "PC-9") # egfr mut


stopifnot(gsub("\\.bed$", "", basename(basename(all_files))) %in% names(bed2annot))

filename = "encode_data/ENCFF001WUY.bed"

if(buildTable) {
    
    all_histOverlap_dt <- foreach(filename=all_files, .combine='rbind') %dopar% {
        
        ds_name <- gsub("\\.bed$", "", basename(basename(filename)))
        stopifnot(ds_name %in% names(bed2annot))
        
        hist_dt <- read.delim(filename, header=F, stringsAsFactors = FALSE,
                              col.names=c("chromo", "chromStart", "chromEnd", "name",
                                          "score","strand","signalValue", "pValue_log10", "qValue_log10", "peak"))
        #chr1    564540  564690  .       0       .       7       5.47141 -1      -1
        #chr1    569820  569970  .       0       .       10      5.00553 -1      -1
        hist_dt$signalValue_qtNorm <- hist_dt$signalValue/quantile(hist_dt$signalValue, probs = normQt)
        hist_dt$signalValue_zNorm <- as.numeric(scale(hist_dt$signalValue, center = TRUE, scale = TRUE))
        
        
        hicds_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
            
            cat(paste0("... start\t", ds_name, "\t", hicds, "\n"))
            
            g2t_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt")
            g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
            g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
            tad_g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
            stopifnot(!duplicated(g2t_dt$entrezID))
            entrez2tad <- setNames(g2t_dt$region, g2t_dt$entrezID)
            
            stopifnot(entrez_oi %in% names(entrez2tad))
            
            tad_oi <- entrez2tad[names(entrez2tad) == entrez_oi]
            
            
            tad_dt <- read.delim(file.path("..", hicds, "genes2tad/all_assigned_regions.txt"), 
                                 col.names=c("chromo", "region", "start", "end"), header=F, stringsAsFactors=FALSE)
            
            stopifnot(tad_oi %in% tad_dt$region)
            tadChr_oi <- tad_dt$chromo[tad_dt$region == tad_oi]
            tadStart_oi <- tad_dt$start[tad_dt$region == tad_oi]
            tadEnd_oi <- tad_dt$end[tad_dt$region == tad_oi]
            
            chromo_dt <- hist_dt[hist_dt$chromo == tadChr_oi,]
            
            
            # overlap: if start or end of the peak within the tad region
            overlap_dt <- chromo_dt[ (chromo_dt$chromStart >= tadStart_oi & chromo_dt$chromStart <= tadEnd_oi) | 
                                         (chromo_dt$chromEnd >= tadStart_oi & chromo_dt$chromEnd <= tadEnd_oi) ,]
            
            
            
            if(all(overlap_dt$pValue_log10) == -1) {
                nSignifPvals <- NA
            } else {
                # if was 3 -> should retrieve 0.001
                nSignifPvals <- sum(10^-(overlap_dt$pValue_log10) <= pvalThresh)
            }
            
            
            if(all(overlap_dt$qValue_log10) == -1) {
                nSignifQvals <- NA
            } else {
                nSignifQvals <- sum(10^-(overlap_dt$qValue_log10) <= qvalThresh)
            }
            
            mean_signal_raw <- mean(overlap_dt$signalValue)
            mean_signal_qt <- mean(overlap_dt$signalValue_qtNorm)
            mean_signal_z <- mean(overlap_dt$signalValue_zNorm)
            
            data.frame(
                hicds = hicds,
                gene_oi=gene_oi,
                tad_oi = tad_oi,
                tadStart_oi=tadStart_oi,
                tadEnd_oi=tadEnd_oi,
                hist_id = ds_name,
                hist_mark = bed2annot[ds_name],
                hist_cl = bed2cl[ds_name],
                mean_signal_raw = mean_signal_raw,
                mean_signal_qt = mean_signal_qt,
                mean_signal_z = mean_signal_z,
                nSignifPvals=nSignifPvals,
                nSignifQvals=nSignifQvals,
                stringsAsFactors = FALSE
            )
            
        }
        hicds_dt
    }
    outFile <- file.path(outFolder, paste0(gene_oi, "_all_histOverlap_dt.Rdata"))
    save(all_histOverlap_dt, file=outFile, version=2)
    cat(paste0("... written: ", outFile, "\n"))
} else {
    # all_histOverlap_dt <- get(load("METH_DATA_FOR_TAD/AKR1C1_all_histOverlap_dt.Rdata"))
    outFile <- file.path(outFolder, paste0(gene_oi, "_all_histOverlap_dt.Rdata"))
    all_histOverlap_dt <- get(load(outFile))
}


all_hists <- unique(all_histOverlap_dt$hist_mark)
hist_m=all_hists[1]

plot_vars <- c("mean_signal_raw", "mean_signal_qt", "mean_signal_z", "nSignifPvals", "nSignifQvals")
plot_var = plot_vars[1]

for(hist_m in all_hists) {
    
    sub_dt <- all_histOverlap_dt[all_histOverlap_dt$hist_mark == hist_m,]

    for(plot_var in plot_vars) {
        
        plotTit <- paste0(hist_m, " - ", plot_var)
        mySub <- paste0("")
        legTitle <- ""
        
        sub_dt$hist_cl <- factor(sub_dt$hist_cl, levels=cl_levels)
        stopifnot(!is.na(sub_dt$hist_cl))
        

        p3 <- ggboxplot(sub_dt,
                        x = "hist_cl",
                        add = "jitter",
                        y = paste0(plot_var),
                        # combine = TRUE,                  # Combine the 3 plots
                        xlab ="",
                        ylab = paste0(plot_var),
                        # add = "median",                  # Add median line.
                        rug = FALSE,                      # Add marginal rug
                        color = "hicds",
                        # fill = "hicds",
                        palette = "jco"
        ) +
            ggtitle(plotTit, subtitle = mySub)+
            # scale_color_manual(values=my_cols)+
            # scale_fill_manual(values=my_cols)  +
            labs(color=paste0(legTitle),fill=paste0(legTitle)) +
            # guides(color=FALSE)+
            scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
            # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
            theme(
                # text = element_text(family=fontFamily),
                panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
                panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
                panel.background = element_rect(fill = "transparent"),
                panel.grid.major.x =  element_blank(),
                panel.grid.minor.x =  element_blank(),
                axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
                axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
                axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
                axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
                # axis.text.x = element_blank(),
                plot.title = element_text(hjust=0.5, size = 16, face="bold"),
                plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
                legend.title = element_text(face="bold")
            ) 
        
        outFile <- file.path(outFolder, paste0(hist_m, "_", plot_var, "_allDS_boxplot.", plotType))
        ggsave(p3, file=outFile, height=myHeight, width=myWidth)
        cat(paste0("... written: ", outFile, "\n"))
        
        
        
        
    }
        
    
}


# Rscript meth_data_for_tad_bw.R 

require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)


plotType <-  "svg"
myHeight <- 7
myWidth <- 10
plotCex <- 1.4


outFolder <- file.path("METH_DATA_FOR_TAD_BW")
dir.create(outFolder, recursive = TRUE)

encode_annot_dt <- read.delim("metadata_plot_list.csv", header=FALSE, sep=",", stringsAsFactors = FALSE)
bed2annot <- setNames(as.character(encode_annot_dt$V3), as.character(encode_annot_dt$V1))
bed2annot <- gsub("-human", "",  bed2annot)
bed2cl <- setNames(as.character(encode_annot_dt$V2), as.character(encode_annot_dt$V1))

buildTable <- TRUE

cl_levels <- c("lung",
               "upper lobe of left lung",
               "AG04450", # fetal
               "IMR-90" , # normal
               "PC-9") # egfr mut


bwOverTAD_folder <- file.path("bwOverTAD")
all_overTAD_files <- list.files(bwOverTAD_folder,pattern=".tab$")
i_file=all_overTAD_files[1]


if(buildTable){
    
    all_histOverTAD_dt <- foreach(i_file=all_overTAD_files, .combine='rbind') %dopar% {
        
        ds_name <- gsub("_overTAD.tab$", "", i_file)
        stopifnot(ds_name %in% names(bed2annot))
        stopifnot(ds_name %in% names(bed2cl))
        
        i_dt <- read.delim(file.path(bwOverTAD_folder,i_file),  header=F, stringsAsFactors = F,
                           col.names=c("name", "size", "covered", "sum", "mean0","mean"))
        i_dt$hist_id <- ds_name
        i_dt$hist_mark <- bed2annot[ds_name]
        i_dt$hist_cl <- bed2cl[ds_name]
        i_dt
    }
    
    
    
    
    outFile <- file.path(outFolder, paste0("all_histOverTAD_dt.Rdata"))
    save(all_histOverTAD_dt, file=outFile, version=2)
    cat(paste0("... written: ", outFile, "\n"))
} else {
    # all_histOverTAD_dt <- get(load("METH_DATA_FOR_TAD_BW_G38/all_histOverTAD_dt.Rdata"))
    outFile <- file.path(outFolder, paste0("all_histOverTAD_dt.Rdata"))
    all_histOverTAD_dt <- get(load(outFile))
}


all_hists <- c("H3K27me3", "H3K27ac")
hist_m <- all_hists[]
plot_vars <- c("mean", "mean0", "sum")
plot_var = plot_vars[1]

for(hist_m in all_hists) {

    sub_dt <- all_histOverTAD_dt[all_histOverTAD_dt$hist_mark == hist_m,]

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
                        # color = "hicds",
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


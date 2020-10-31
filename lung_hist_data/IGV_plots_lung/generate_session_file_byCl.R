require(ggsci)
# Rscript ./generate_session_file_byCl.R

files_dt <- read.delim("../metadata_plot_list.csv", stringsAsFactors = FALSE, sep=",",
                       col.names=c("data", "origin", "histmark"), header=F)
files_dt <- files_dt[order(files_dt$histmark, files_dt$origin),]

topTADfile <- "../AKR1C_TAD_4880001_5120000.bed"

data_folder <- "../encode_data"

plot_start <- paste0(4880000 - 10000)#"4870000"  # 4880000	5120000	
plot_end <- paste0(5150000 + 10000)#"51320001"
# True TAD: 


cl_types <- unique(files_dt$origin)
cl_types <- c( "PC-9"  , "IMR-90", "lung"   ,  "upper lobe of left lung" ,"AG04450" )
stopifnot(setequal(cl_types, unique(files_dt$origin)))

cl_colors <- setNames(pal_jama()(length(cl_types)), cl_types)
cl_colors_rgb <- sapply(cl_colors, function(x)paste0(as.vector(col2rgb(x)), collapse=","))

# [1] "upper lobe of left lung" "IMR-90"                  "lung"                    "PC-9"                   
# [5] "AG04450"                


#histMark="H3K27ac"
#histMark="H3K4me1"

#bwFolder="/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017"

#hicds="GSE105318_DLD1_40kb"
#hicds="ENCSR504OTV_transverse_colon_40kb"
#hicds="Rao_HCT-116_2017_40kb"
#topTADfile="/media/electron/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/IGV_PLOTS/${hicds}_TCGAcoad_msi_mss_TADcoord.bed"

# refPeaksFile="/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE77737_Cohen2017/$histMark/${histMark}_all_merged_peaks_cutAndSorted_filtered_merged_named.bed"

files_dt_init <- files_dt

for(cl in cl_types) {
  
  files_dt <- files_dt_init[files_dt_init$origin == cl,]
  stopifnot(nrow(files_dt) > 0)
 
  outFile <- file.path(paste0(gsub(" ", "", cl), "_lung_data_IGVsession.xml"))
  file.remove(outFile)
  
  
  fileConn <- file(outFile, "w")
  writeLines(text="<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>
  <Session genome=\"hg19\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" path=\"/home/marie/Desktop/29.06/IGV_session_all_samples/GSE77737_MSI_MSS_igv_session_H3K427ac.xml\" version=\"8\">
  <Resources>",con=fileConn)
  writeLines(text=paste0("\t<Resource path=\"", topTADfile, "\"/>"), con=fileConn)
  writeLines(text=paste0("\t<Resource name=\"Refseq Genes\" path=\"https://s3.dualstack.us-east-1.amazonaws.com/igv.org.genomes/hg19/ncbiRefGene.txt.gz\"/>"), con=fileConn)	
  
  all_files <- file.path(data_folder, paste0(files_dt$data, ".bed"))
  
  for(histf in all_files) {
    writeLines(text=paste0("\t<Resource path=\"", histf, "\"/>"), con=fileConn)
  }
  
  writeLines(text="</Resources>
  <Panel height=\"361\" name=\"DataPanel\" width=\"2476\">", con=fileConn)
  
  
  
  for(i in 1:nrow(files_dt)) {
    track_name <- paste0(files_dt$origin[i], "_", files_dt$histmark[i],  "_", files_dt$data[i])
    track_color <- cl_colors_rgb[files_dt$origin[i]]
    track_path <- file.path(data_folder, paste0(files_dt$data[i], ".bed"))
    
    cat(paste0(track_path,"\n"))
    
    writeLines(text=paste0("\t\t<Track altColor=\"0,0,178\" autoScale=\"false\" clazz=\"org.broad.igv.track.DataSourceTrack\" color=\"", track_color,
                           "\" displayMode=\"COLLAPSED\" featureVisibilityWindow=\"-1\" fontSize=\"10\" id=\"", track_path, 
                           "\" name=\"", track_name, "\" normalize=\"false\" renderer=\"BAR_CHART\" sortable=\"true\" visible=\"true\" windowFunction=\"mean\">"), con=fileConn)
    
    writeLines(text= paste0( "\t\t\t<DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"3.1860125\" minimum=\"0.0\" type=\"LINEAR\"/>"), con=fileConn)
    
    writeLines(text="\t\t</Track>", con=fileConn)
    
    
  }
  
  writeLines(text = paste0("</Panel> 
  <Panel height=\"837\" name=\"FeaturePanel\" width=\"2476\">
        <Track attributeKey=\"Reference sequence\" clazz=\"org.broad.igv.track.SequenceTrack\" fontSize=\"10\" id=\"Reference sequence\" name=\"Reference sequence\" visible=\"true\"/>
        <Track attributeKey=\"Refseq Genes\" clazz=\"org.broad.igv.track.FeatureTrack\" color=\"0,0,178\" colorScale=\"ContinuousColorScale;0.0;444.0;255,255,255;0,0,178\" fontSize=\"10\" id=\"https://s3.dualstack.us-east-1.amazonaws.com/igv.org.genomes/hg19/ncbiRefGene.txt.gz\" name=\"Refseq Genes\" visible=\"true\"/>
  <Track altColor=\"0,0,178\" autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" color=\"0,0,178\" displayMode=\"COLLAPSED\" featureVisibilityWindow=\"-1\" fontSize=\"10\" id=\"", topTADfile, "\" name=\"topTADs\" renderer=\"BASIC_FEATURE\" sortable=\"false\" visible=\"true\" windowFunction=\"count\"/>
  </Panel>
  <PanelLayout dividerFractions=\"0.3012448132780083\"/>
  <HiddenAttributes>
  <Attribute name=\"DATA FILE\"/>
  <Attribute name=\"DATA TYPE\"/>
  <Attribute name=\"NAME\"/>
  </HiddenAttributes>
  </Session>"), con=fileConn)
  
  
  
  cat(paste0("WRITTEN: ", outFile, "\n"))
  
  outFile2 <- file.path(paste0(gsub(" ", "", cl), "_lung.batch"))
  file.remove(outFile2)
  
  fileConn2 <- file(outFile2, "w")
  
  writeLines(
  paste0("snapshotDirectory /media/electron/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/lung_hist_data/IGV_plots_lung/plots\n
  maxPanelHeight 1200\n
  goto chr10:", plot_start, "-", plot_end,"\n",
  "snapshot ", gsub(" ", "_", cl), "_AKR1C_TAD_chr10_", plot_start, "_", plot_end, "_slop1000.png"), con=fileConn2)
  
  cat(paste0("WRITTEN: ", outFile2, "\n"))
  
  
  
  
   
}


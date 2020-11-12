require(ggsci)
# Rscript ./generate_session_file.R

# outSetDir 

setDir <- "/media/electron"
# setDir <- ""
dataFolder <- file.path(setDir, "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/daniele_lung_hist_data/SUBSET_FOR_IGV")
filePatt <- "withImr90.+FinalNorm.bedgraph"

plotStart <- "4835000"
plotEnd <- "5080000"

plotStart="4820000"
plotEnd="5108000"

plotStart="4800000"
plotEnd="5110000"

filePrefix <- paste0("subset_chr10_", plotStart, "_", plotEnd, "_marie_withImr90")
fileSuffix <- "FinalNorm.bedgraph"

topTADfile <- "AKR1C_TAD_4837809_5077808.bed"

all_files <- list.files(dataFolder, full.names = TRUE, pattern=paste0(filePrefix, ".+", fileSuffix))

cl_types <- gsub(paste0(filePrefix, "_(.+)_", fileSuffix), "\\1", basename(all_files))


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

outFile <- file.path(paste0("lung_data_IGVsession_", plotStart, "_", plotEnd, ".xml"))
file.remove(outFile)

# msi_color="255,0,0"
# mss_color="0,178,0"
# refpeak_color="255,140,0"
# 
# msi_bwFiles=($(ls $bwFolder/MSI/$histMark/*_${histMark}.bw ))
# mss_bwFiles=($(ls $bwFolder/MSS/$histMark/*_${histMark}.bw ))


# for bw_file in ${msi_bwFiles[@]}; do
# 	echo $bw_file
# done

# cat > $outFile <<- EOM
# 	<?xml version="1.0" encoding="UTF-8" standalone="no"?>
# 	<Session genome="hg19" hasGeneTrack="true" hasSequenceTrack="true" path="/home/marie/Desktop/29.06/IGV_session_all_samples/GSE77737_MSI_MSS_igv_session_H3K427ac.xml" version="8">
# 		<Resources>
# EOM

fileConn <- file(outFile, "w")
writeLines(text="<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>
  <Session genome=\"hg38\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" path=\"/home/marie/Desktop/29.06/IGV_session_all_samples/GSE77737_MSI_MSS_igv_session_H3K427ac.xml\" version=\"8\">
  <Resources>",con=fileConn)
writeLines(text=paste0("\t<Resource path=\"", topTADfile, "\"/>"), con=fileConn)
writeLines(text=paste0("\t<Resource name=\"Refseq Genes\" path=\"https://s3.dualstack.us-east-1.amazonaws.com/igv.org.genomes/hg38/ncbiRefGene.txt.gz\"/>"), con=fileConn)	

for(histf in all_files) {
  writeLines(text=paste0("\t<Resource path=\"", histf, "\"/>"), con=fileConn)
}

writeLines(text="</Resources>
  <Panel height=\"361\" name=\"DataPanel\" width=\"2476\">", con=fileConn)
  

# for bw_file in ${msi_bwFiles[@]}; do
#         echo -e "\t<Resource path=\"$bw_file\"/>" >> $outFile
# done
# echo -e "\t<Resource path=\"$topTADfile\"/>" >> $outFile
# echo -e "\t<Resource path=\"$refPeaksFile\"/>" >> $outFile
# for bw_file in ${mss_bwFiles[@]}; do
#         echo -e "\t<Resource path=\"$bw_file\"/>" >> $outFile
# done
# 
# cat >> $outFile <<- EOM
#     </Resources>
#     <Panel height="361" name="DataPanel" width="2476">
# EOM

# for(i in 1:nrow(files_dt)) {

  for(myfile in all_files) {
  
  track_name <- gsub("marie_", "", basename(myfile)) 
  track_name <- gsub("\\.bedgraph", "", track_name)
  track_color <- cl_colors_rgb[gsub(paste0(filePrefix, "_(.+)_", fileSuffix), "\\1", basename(myfile))]
  stopifnot(!is.na(track_color))
  track_path <- myfile
  
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
        <Track attributeKey=\"Refseq Genes\" clazz=\"org.broad.igv.track.FeatureTrack\" color=\"0,0,178\" colorScale=\"ContinuousColorScale;0.0;444.0;255,255,255;0,0,178\" fontSize=\"10\" id=\"https://s3.dualstack.us-east-1.amazonaws.com/igv.org.genomes/hg38/ncbiRefGene.txt.gz\" name=\"Refseq Genes\" visible=\"true\"/>
  <Track altColor=\"0,0,178\" autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" color=\"0,0,178\" displayMode=\"COLLAPSED\" featureVisibilityWindow=\"-1\" fontSize=\"10\" id=\"", topTADfile, "\" name=\"topTADs\" renderer=\"BASIC_FEATURE\" sortable=\"false\" visible=\"true\" windowFunction=\"count\"/>
  </Panel>
  <PanelLayout dividerFractions=\"0.3012448132780083\"/>
  <HiddenAttributes>
  <Attribute name=\"DATA FILE\"/>
  <Attribute name=\"DATA TYPE\"/>
  <Attribute name=\"NAME\"/>
  </HiddenAttributes>
  </Session>"), con=fileConn)

# <Track altColor=\"$refpeak_color\" autoScale=\"false\" clazz=\"org.broad.igv.track.FeatureTrack\" color=\"$refpeak_color\" displayMode=\"COLLAPSED\" featureVisibilityWindow=\"-1\" fontSize=\"10\" id=\"$refPeaksFile\" name=\"refPeaks\" renderer=\"BASIC_FEATURE\" sortable=\"false\" visible=\"true\" windowFunction=\"count\"/>

  
cat(paste0("WRITTEN: ", outFile, "\n"))


# for bw_file in ${msi_bwFiles[@]}; do
# 	tmp=`basename $bw_file`
# 	track_name="${histMark}_$tmp"
# 	echo -e "\t\t<Track altColor=\"0,0,178\" autoScale=\"false\" clazz=\"org.broad.igv.track.DataSourceTrack\" color=\"$msi_color\" displayMode=\"COLLAPSED\" featureVisibilityWindow=\"-1\" fontSize=\"10\" id=\"$bw_file\" name=\"$track_name\" normalize=\"false\" renderer=\"BAR_CHART\" sortable=\"true\" visible=\"true\" windowFunction=\"mean\">" >> $outFile
#     echo -e "\t\t\t<DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"3.1860125\" minimum=\"0.0\" type=\"LINEAR\"/>" >> $outFile
#     echo -e "\t\t</Track>" >> $outFile
# done
# for bw_file in ${mss_bwFiles[@]}; do
# 	tmp=`basename $bw_file`
# 	track_name="${histMark}_$tmp"
# 	echo -e "\t\t<Track altColor=\"0,0,178\" autoScale=\"false\" clazz=\"org.broad.igv.track.DataSourceTrack\" color=\"$mss_color\" displayMode=\"COLLAPSED\" featureVisibilityWindow=\"-1\" fontSize=\"10\" id=\"$bw_file\" name=\"$track_name\" normalize=\"false\" renderer=\"BAR_CHART\" sortable=\"true\" visible=\"true\" windowFunction=\"mean\">" >> $outFile
#     echo -e "\t\t\t<DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"3.1860125\" minimum=\"0.0\" type=\"LINEAR\"/>" >> $outFile
#     echo -e "\t\t</Track>" >> $outFile
# done

# cat >> $outFile <<- EOM
#     </Panel>
#     <Panel height="837" name="FeaturePanel" width="2476">
#         <Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" sortable="false" visible="true"/>
#         <Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;423.0;255,255,255;0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="35" id="hg19_genes" name="RefSeq Genes" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count">
#             <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="423.0" minimum="0.0" type="LINEAR"/>
#         </Track>
#         <Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="$topTADfile" name="topTADs" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>
#         <Track altColor="$refpeak_color" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="$refpeak_color" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="$refPeaksFile" name="refPeaks" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>
#     </Panel>
#     <PanelLayout dividerFractions="0.3012448132780083"/>
#     <HiddenAttributes>
#         <Attribute name="DATA FILE"/>
#         <Attribute name="DATA TYPE"/>
#         <Attribute name="NAME"/>
#     </HiddenAttributes>
# </Session>
# 
# EOM
# 
# echo "WRITTEN: $outFile"

# 
# ###################################################################################################################################################
# ########## END ####################################################################################################################################
# ###################################################################################################################################################
# echo "*** DONE"
# echo "... written: $logFile"
# end_time=$(date -R)    
# echo $start_time
# echo $end_time
# exit 0


http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/7.0/c3.tft.v7.0.entrez.gmt
http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/7.0/c3.mir.v7.0.entrez.gmt
http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/7.0/c3.all.v7.0.entrez.gmt

=> dld from http://software.broadinstitute.org/gsea/msigdb/collections.jsp, 11.12.2019


http://www.tfcheckpoint.org/data/TFCheckpoint_download_180515.txt , 11.12.2019

https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv, 11.12.2019

https://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/encodetfppi/gene_set_library_crisp.gmt.gz


http://tfbsdb.systemsbiology.net/download 12.12.2019


HUMAN_hg19_BBLS_1_00_FDR_0_10.bed dld from http://www.igb.uci.edu/~motifmap/motifmap/HUMAN/hg19/multiz46way_placental/HUMAN.hg19_multiz46way.tar.bz2 13.12.2019

https://www.biostars.org/p/164617/
http://rest.kegg.jp/link/pathway/hsa, mv hsa hsa_kegg_entrez.txt dld 13.12.2019


https://amp.pharm.mssm.edu/chea3/ # dld 31.01.2020
wget https://amp.pharm.mssm.edu/chea3/assets/tflibs/lung.TFs.gmt
wget https://amp.pharm.mssm.edu/chea3/assets/tflibs/all_tissues.TFs.gmt



## enhanceratlas data - 14.02.2020
wget http://www.enhanceratlas.org/data/download/enhancer/hs/Lung.bed --output-file enhanceratlas_enhancer_Lung.bed
wget http://www.enhanceratlas.org/data/AllEPs/hs/Lung_EP.txt --output-file enhanceratlas_AllEPs_Lung_EP.bed

Data format EP file
Data format in the file listed here (10 columns):
chrom-Enh - Name of the chromosome for enhancer.
chromStart - The starting position of enhancer.
chromEnd - The ending position of enhancer.
Gene ID - A number taken as the gene id.
chrom-Gene - Name of the chromosome for gene.
TSS - the position of the gene transcription start site.
Transcript - The gene Transcript
NULL
signalValue - Measurement of average enrichment for enhancer.
EP Score - The confidence score of the enhancer-target interaction.

http://www.enhanceratlas.org/Data_format_EP.txt

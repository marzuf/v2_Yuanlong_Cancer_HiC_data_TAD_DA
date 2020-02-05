grep  hg18 inferCARs_data/Orthology.Blocks > inferCARs_data/Orthology.Blocks_hg18.txt

hg18.chr1:3522431-3987024 + [2] (0)


sed -n "s/hg18.\(chr.\+\):\(.\+\)-\(.\+\) +.\+/\1:\2-\3/p" inferCARs_data/Orthology.Blocks_hg18.txt | cut -f1,2,3 > inferCARs_data/Orthology.Blocks_hg18.txt.bed

https://genome.ucsc.edu/cgi-bin/hgLiftOver => converted hg18 to hg19 the file Orthology.Blocks_hg18.txt.bed [with default parameters] -> Orthology.Blocks_hg18.txt.bed_hg19conv.bed

failure file: convert_failure.txt -> save the bed that have not been converted

grep chr inferCARs_data/convert_failure.txt > inferCARs_data/convert_failure_bed.txt

perl process_inferCARs_orthologyBlocks.pl


Rscript convert_hg19_processed_inferCARs_orthology_blocks.R

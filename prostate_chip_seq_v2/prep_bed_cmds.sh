cp ../GSE118514_22Rv1_40kb/genes2tad/all_assigned_regions.txt GSE118514_22Rv1_40kb_all_assigned_
regions.txt

cp ../GSE118514_RWPE1_40kb/genes2tad/all_assigned_regions.txt GSE118514_RWPE1_40kb_all_assigned_
regions.txt
22Rv1_ENCFF282RLY_fcOverControl_rep12.bigWig
22Rv1_ENCFF286BKT_fcOverControl_rep2.bigWig
22Rv1_ENCFF408PLR_fcOverControl_rep1.bigWig
bigwigCompare_22Rv1_RWPE1.bigWig
RWPE1_ENCFF039XYU_fcOverControl_rep12.bigWig
RWPE1_ENCFF053MIK_fcOverControl_rep2.bigWig
RWPE1_ENCFF803NYF_fcOverControl_rep1.bigWig

/mnt/ed4/marie/software/bigWigInfo RWPE1_ENCFF053MIK_fcOverControl_rep2.bigWig
/mnt/ed4/marie/software/bigWigInfo 22Rv1_ENCFF408PLR_fcOverControl_rep1.bigWig

marie@electron:/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/prostate_chip_seq_v2$ 


awk '{ print $1"\t"$3"\t"$4"\t"$2 }' ../GSE118514_RWPE1_40kb/genes2tad/GSE118514_RWPE1_40kb_all_assigned_regions.txt > GSE118514_RWPE1_40kb_all_assigned_regions_ordrd.txt

awk '{ print $1"\t"$3"\t"$4"\t"$2 }' ../GSE118514_22Rv1_40kb/genes2tad/GSE118514_22Rv1_40kb_all_assigned_regions.txt > GSE118514_22Rv1_40kb_all_assigned_regions_ordrd.txt

########### convert bigwig to bedgraph
 /mnt/ed4/marie/software/bigWigToBedGraph ../prostate_chip_seq/22Rv1_ENCFF282RLY_fcOverControl_rep12.bigWig ../prostate_chip_seq/22Rv1_ENCFF282RLY_fcOverControl_rep12.bedGraph
 /mnt/ed4/marie/software/bigWigToBedGraph ../prostate_chip_seq/22Rv1_ENCFF286BKT_fcOverControl_rep2.bigWig ../prostate_chip_seq/22Rv1_ENCFF286BKT_fcOverControl_rep2.bedGraph
 /mnt/ed4/marie/software/bigWigToBedGraph ../prostate_chip_seq/22Rv1_ENCFF408PLR_fcOverControl_rep1.bigWig ../prostate_chip_seq/22Rv1_ENCFF408PLR_fcOverControl_rep1.bedGraph
 /mnt/ed4/marie/software/bigWigToBedGraph ../prostate_chip_seq/RWPE1_ENCFF039XYU_fcOverControl_rep12.bigWig ../prostate_chip_seq/RWPE1_ENCFF039XYU_fcOverControl_rep12.bedGraph
 /mnt/ed4/marie/software/bigWigToBedGraph  ../prostate_chip_seq/RWPE1_ENCFF053MIK_fcOverControl_rep2.bigWig ../prostate_chip_seq/RWPE1_ENCFF053MIK_fcOverControl_rep2.bedGraph
 /mnt/ed4/marie/software/bigWigToBedGraph ../prostate_chip_seq/RWPE1_ENCFF803NYF_fcOverControl_rep1.bigWig ../prostate_chip_seq/RWPE1_ENCFF803NYF_fcOverControl_rep1.bedGraph
 

########### average bigwig over selected TADs
/mnt/ed4/marie/software/bigWigAverageOverBed ../prostate_chip_seq/RWPE1_ENCFF039XYU_fcOverControl_rep12.bigWig ../prostate_chip_seq/GSE118514_RWPE1_40kb_all_assigned_regions_ordrd.txt chip_RWPE1_cover_RWPE1_TADs_rep12.bed
/mnt/ed4/marie/software/bigWigAverageOverBed ../prostate_chip_seq/RWPE1_ENCFF803NYF_fcOverControl_rep1.bigWig ../prostate_chip_seq/GSE118514_RWPE1_40kb_all_assigned_regions_ordrd.txt chip_RWPE1_cover_RWPE1_TADs_rep1.bed
/mnt/ed4/marie/software/bigWigAverageOverBed ../prostate_chip_seq/RWPE1_ENCFF053MIK_fcOverControl_rep2.bigWig ../prostate_chip_seq/GSE118514_RWPE1_40kb_all_assigned_regions_ordrd.txt chip_RWPE1_cover_RWPE1_TADs_rep2.bed

/mnt/ed4/marie/software/bigWigAverageOverBed ../prostate_chip_seq/22Rv1_ENCFF282RLY_fcOverControl_rep12.bigWig ../prostate_chip_seq/GSE118514_22Rv1_40kb_all_assigned_regions_ordrd.txt chip_22Rv1_cover_22Rv1_TADs_rep12.bed
/mnt/ed4/marie/software/bigWigAverageOverBed ../prostate_chip_seq/22Rv1_ENCFF408PLR_fcOverControl_rep1.bigWig ../prostate_chip_seq/GSE118514_22Rv1_40kb_all_assigned_regions_ordrd.txt chip_22Rv1_cover_22Rv1_TADs_rep1.bed
/mnt/ed4/marie/software/bigWigAverageOverBed ../prostate_chip_seq/22Rv1_ENCFF286BKT_fcOverControl_rep2.bigWig ../prostate_chip_seq/GSE118514_22Rv1_40kb_all_assigned_regions_ordrd.txt chip_22Rv1_cover_22Rv1_TADs_rep2.bed



/mnt/ed4/marie/software/bigWigAverageOverBed ../prostate_chip_seq/RWPE1_ENCFF039XYU_fcOverControl_rep12.bigWig ../prostate_chip_seq/GSE118514_22Rv1_40kb_all_assigned_regions_ordrd.txt chip_RWPE1rep12_cover_22Rv1_TADs.bed


/mnt/ed4/marie/software/bigWigAverageOverBed ../prostate_chip_seq/22Rv1_ENCFF282RLY_fcOverControl_rep12.bigWig ../prostate_chip_seq/GSE118514_RWPE1_40kb_all_assigned_regions_ordrd.txt chip_22Rv1rep12_cover_RWPE1_TADs.bed


# /mnt/ed4/marie/software/bigWigAverageOverBed tmp.bw tmp_tad.bed tmp_cover.tad
#output is:
#   name - name field from bed, which should be unique
#   size - size of bed (sum of exon sizes
#   covered - # bases within exons covered by bigWig
#   sum - sum of values over all bases covered
#   mean0 - average over bases with non-covered bases counting as zeroes
#   mean - average over just covered bases
  



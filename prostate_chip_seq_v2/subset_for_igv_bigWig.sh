
# ./subset_for_igv_bigWig.sh

# 22Rv1_ENCFF286BKT_fcOverControl.bedGraph
# RWPE1_ENCFF039XYU_fcOverControl.bedGraph

# top3 rwpe1 tads
#chr12	54160001	54440000	chr12_TAD194
#chr7	116080001	116320000	chr7_TAD424
#chr17	46720001	46880000	chr17_TAD174

dataFolder="../prostate_chip_seq"

filePrefix=""
fileSuffix="_fcOverControl.bedGraph"

# 4837809_5077808.
#tadChromo="chr10"
##tadStart="4820000"
##tadEnd="5108000"
#tadStart="4837809"
#tadEnd="5077808"
##tadStart="0"
##tadEnd="200"

all_chromos=( "chr12" "chr7" "chr17" )
all_starts=( 54160001 116080001 46720001 )
all_ends=( 54440000 116320000 46880000 )

outFolder="SUBSET_FOR_IGV_BIGWIG"
mkdir -p $outFolder


for i in "${!all_chromos[@]}"; do
	
	
	tadChromo=${all_chromos[$i]}	
	tadStart=${all_starts[$i]}	
	tadEnd=${all_ends[$i]}	
	
	echo "$tadChromo - $tadStart - $tadEnd"		
	 
	for f in 22Rv1_ENCFF282RLY_rep12 22Rv1_ENCFF286BKT_rep2 22Rv1_ENCFF408PLR_rep1 RWPE1_ENCFF039XYU_rep12 RWPE1_ENCFF053MIK_rep2 RWPE1_ENCFF803NYF_rep1; do

		infile="$dataFolder/${filePrefix}${f}${fileSuffix}"
	#	infile="tmp.bedgraph"

		outfile="$outFolder/subset_${tadChromo}_${tadStart}_${tadEnd}${filePrefix}_${f}${fileSuffix}"

		awk -v chrom=$tadChromo -v start=$tadStart -v end=$tadEnd '{if($1 == chrom && $2 >= start && $3 <= end ) {print}}' ${infile} > $outfile

		#echo awk -v chrom=$tadChromo -v start=$tadStart -v end=$tadEnd {if($1 == chrom && $2 >=start && $3 <=end ){print}} ${infile} #> $outfile
		#echo awk '{if($1 =='$tadChromo && $2 >=$tadStart && $3 <=$tadEnd ){print}}' ${infilef}
		# > "$outFolder/subset_${tadChromo}_${tadStart}_${tadEnd}_/${filePrefix}_${f}_${fileSuffix}"
		
		echo $outfile
	done
done

#!/usr/bin/bash

all_files=( $( ls *_40kb/TCGA*/*.txt ) )
#my_array=( $( my_command) )

#    Column1: Unique Peak ID
#    Column2: chromosome
#    Column3: starting position
#    Column4: ending position
#    Column5: Strand (+/- or 0/1, where 0="+", 1="-")



#chr10   chr10_BOUND1    1       80000
#chr10   chr10_TAD1      80001   480000
#chr10   chr10_TAD2      480001  840000

homer_script="bin/findMotifsGenome.pl"
homer_genome="hg19r" # r -> --mask
#findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
# i.e. findMotifsGenome.pl ERpeaks.txt hg18 ER_MotifOutput/ -size 200 -mask 

ncpu=80

#PATH="bin"


#HMEC_40kb/TCGAbrca_lum_bas/HMEC_40kb_TCGAbrca_lum_bas_adjPvalComb_0.01_minus.txt
#HMEC_40kb/TCGAbrca_lum_bas/HMEC_40kb_TCGAbrca_lum_bas_adjPvalComb_0.01_plus.txt
#HMEC_40kb/TCGAbrca_lum_bas/HMEC_40kb_TCGAbrca_lum_bas_signifFDR_0.2_minus.txt
#HMEC_40kb/TCGAbrca_lum_bas/HMEC_40kb_TCGAbrca_lum_bas_signifFDR_0.2_plus.txt

allsignifs=( "adjPvalComb_0.01" "signifFDR_0.2")

for file in ${all_files[@]}; do

#    grep _TAD $file | cut -f2,1,3,4,2 > ${file}_plus.txt

    x1=`dirname $file`
    hicds=`dirname $x1`
    exprds=`basename $x1`

    echo "$hicds - $exprds"



    for signif in ${allsignifs[@]}; do




        filePlus="$hicds/$exprds/${hicds}_${exprds}_${signif}_plus.txt"
        fileMinus="$hicds/$exprds/${hicds}_${exprds}_${signif}_minus.txt"
        
        filePlusMinus="$hicds/$exprds/${hicds}_${exprds}_${signif}_plusminus.txt"
        cat $filePlus $fileMinus > $filePlusMinus


        echo "... written: $filePlusMinus"

        mkdir -p $hicds/$exprds/$signif/MotifOutput_plus
        mkdir -p $hicds/$exprds/$signif/MotifOutput_minus
        mkdir -p $hicds/$exprds/$signif/MotifOutput_plusminus

        echo $homer_script $filePlus hg19r $hicds/$exprds/$signif/MotifOutput_plus -size 200 -p $ncpu
        $homer_script $filePlus hg19r $hicds/$exprds/$signif/MotifOutput_plus -size 200 -p $ncpu
        echo $homer_script $fileMinus hg19r $hicds/$exprds/$signif/MotifOutput_minus -size 200 -p $ncpu
        $homer_script $fileMinus hg19r $hicds/$exprds/$signif/MotifOutput_minus -size 200 -p $ncpu
        echo $homer_script $filePlusMinus hg19r $hicds/$exprds/$signif/MotifOutput_plusminus -size 200 -p $ncpu
        $homer_script $filePlusMinus hg19r $hicds/$exprds/$signif/MotifOutput_plusminus -size 200 -p $ncpu

        #exit 0

    done

done








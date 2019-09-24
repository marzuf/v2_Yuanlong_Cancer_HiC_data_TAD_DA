#!/usr/bin/bash

all_files=( $( ls ../*/genes2tad/all_assigned_regions.txt ) )
#my_array=( $( my_command) )


#    Column1: chromosome
#    Column2: starting position
#    Column3: ending position
#    Column4: Unique Peak ID
#    Column5: not used
#    Column6: Strand (+/- or 0/1, where 0="+", 1="-")


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

ncpu=40

#PATH="bin"

for file in ${all_files[@]}; do

#    grep _TAD $file | cut -f2,1,3,4,2 > ${file}_plus.txt

    x1=`dirname $file`
    x2=`dirname $x1`
    x3=`basename $x2`

    #echo $x3
    mkdir $x3

    filePlus="$x3/${x3}_assigned_regions_plus.txt"
    fileMinus="$x3/${x3}_assigned_regions_minus.txt"

    filePlusMinus="$x3/${x3}_assigned_regions_plusminus.txt"


    awk '{print $2"\t"$1"\t"$3"\t"$4"\t+"}' $file | grep _TAD  > $filePlus
    awk '{print $2"\t"$1"\t"$3"\t"$4"\t-"}' $file | grep _TAD  > $fileMinus

    cat $filePlus $fileMinus > $filePlusMinus


    echo "... written: $filePlus"
    echo "... written: $fileMinus"
    echo "... written: $filePlusMinus"

    mkdir $x3/MotifOutput_plus

    $homer_script $filePlus hg19r $x3/MotifOutput_plus -size 200 -p $ncpu
    $homer_script $fileMinus hg19r $x3/MotifOutput_minus -size 200 $ncpu
    $homer_script $filePlusMinus hg19r $x3/MotifOutput_plusminus -size 200


    exit 0

done


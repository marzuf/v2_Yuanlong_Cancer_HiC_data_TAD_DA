# change in the pipeline 13.08.2019

*** TAD extraction from Yuanlong (discussion with Yuanlong, 13.08.2019)
- if TADs have the same CD rank -> they belong to the same TAD if they are not separated by more than 50 bins, e.g.
=> iterate over the rows and check that if CD_rank of row i == CD_rank of row i-1, diff(binIdx of row i, binIdx of row i-1) is <= 50

        binIdx  CD_rank
        15          1
        17          1
        22          2
        ...         ..

- discard also if there is a TAD that contains the mid position of the centromere

*** min count (discussion with Giovanni, 13.08.2019)
- retain a gene if it has at least min_counts counts in at least min_sampleRatio of the samples (-> at least 5 counts in 80% of the samples; in other words, if at least 80% of the samples have 5 or more counts)

# STEPS TO RUN

"0cleanInputTCGAminCount"
"1cleanInput"
"3"
"4"
"5"
"5sameNbr"  
"6" 
"7"
"7sameNbr" 
"8c"         # "8cOnlyRatioDown" # will need all for the wave plot
"9" 
"10"
"10sameNbr"
"11" 
"11sameNbr"
"19sameNbr"
"19onlyFC"
"14f2"


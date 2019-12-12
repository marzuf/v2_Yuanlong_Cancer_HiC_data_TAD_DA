#!/usr/bin/perl -w
use strict;
use warnings;

# perl process_c3_set.pl

my $inFile = "c3.mir.v7.0.entrez.gmt";
my $outFile = "c3.mir.v7.0.entrez_processed.txt";
#my $inFile = "c3.all.v7.0.entrez.gmt";
#my $outFile = "c3.all.v7.0.entrez_processed.txt";
#my $inFile = "c3.tft.v7.0.entrez.gmt";
#my $outFile = "c3.tft.v7.0.entrez_processed.txt";

system("rm -f $outFile");


open(my $outF, ">>", $outFile);

my $header = "regElement\tentrezID";
printf($outF $header."\n");


open IN, $inFile || die;

while(my $line = <IN>) {
    chomp $line;
    my @fields = split(' ', $line);

	my $tf = shift @fields;
	shift @fields; # second field is an URL

    # at the 1st line -> til 1 field
    foreach my $gene (@fields) {
        printf($outF "$tf\t$gene\n");
    }

}



close($outF);

print "!!!!!!!!!! DONE !!!!!!!!!!\n";
print "... written in $outFile\n";



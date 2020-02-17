#!/usr/bin/perl -w
use strict;
use warnings;

# perl process_chea3_encode.pl

my $inFile = "ENCODE_ChIP-seq.gmt";




my $outFile = "chea3_ENCODE_processed.txt";
system("rm -f $outFile");


open(my $outF, ">>", $outFile);

my $header = "regSymbol\ttargetSymbol";
printf($outF $header."\n");


open IN, $inFile || die;

while(my $line = <IN>) {
    chomp $line;
    my @fields = split(' ', $line);

	my $tf = shift @fields;

    # at the 1st line -> til 1 field
    foreach my $gene (@fields) {
        printf($outF "$tf\t$gene\n");
    }
}



close($outF);

print "!!!!!!!!!! DONE !!!!!!!!!!\n";
print "... written in $outFile\n";



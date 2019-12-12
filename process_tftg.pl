#!/usr/bin/perl -w
use strict;
use warnings;

# perl process_tftg.pl

my $inFolder="geneHits_minus5Kbp_plus5Kbp";

my @all_files=`ls $inFolder`;


my $outFile = "tftg_db_all_processed.txt";
system("rm -f $outFile");
open (my $outF, ">>", $outFile);

my $header = "regElement\tentrezID";
printf($outF $header."\n");

foreach my $inFile(@all_files){

	print $inFile."\n";

	# geneHits_minus5Kbp_plus5Kbp/fullGenome_motifHits_Zscan4.2.csv 

	$inFile =~ m/fullGenome_motifHits_(.+).csv/;

	my $tf=$1;

	print $tf."\n"; 

	open IN, "$inFolder/$inFile" || die;

	while(my $line = <IN>) {
		chomp $line;
		my @fields = split(',', $line);
		my $target = shift @fields;
		next if $target eq "Entrez ID";
		printf($outF "$tf\t$target\n");
	}
	#exit;
}

close($outF);
print "!!!!!!!!!! DONE !!!!!!!!!!\n";
print "... written in $outFile\n";




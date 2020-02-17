
use strict;
use warnings;

# perl process_proCARs_orthologyBlocks.pl

my $inputFile = "procars_orthology_blocks.txt";

my $outputFile = "procars_orthology_blocks_processsed.txt";
system("rm -f $outputFile");

#>1
#homo_sapiens.1:2328555-2560923 +
#pan_troglodytes.1:2276241-2516460 +
#pongo_abelii.1:227927143-228157951 -
#macaca_mulatta.1:5465631-5712980 +
#mus_musculus.4:154257282-154455012 -



open(my $outF, ">>", $outputFile);

my $header = "blockID\tgenome\tchromo\tstart\tend";
printf($outF $header."\n");

my $blockID;

open IN, $inputFile || die;

while(my $line = <IN>) {
    chomp $line;

	next if $line =~ /^\s*$/;	# skip white lines

	if ($line =~ /^>/) {

		$blockID = $line;
		next;

	} else {

		$line =~ s/(.+)\.(.+):(\d+)-(\d+) .+/$1\tchr$2\t$3\t$4/g;
		# print $line;
		# exit;
        printf($outF "$blockID\t$line\n");

	}
}



close($outF);

print "!!!!!!!!!! DONE !!!!!!!!!!\n";
print "... written in $outputFile\n";


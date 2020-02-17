
use strict;
use warnings;

# perl process_inferCARs_orthologyBlocks.pl

my $inputFile = "inferCARs_data/Orthology.Blocks";

my $outputFile = "inferCARs_data/Orthology.Blocks_processed.txt";
system("rm -f $outputFile");

#>1
#hg18.chr1:835227-3324958 + [2] (0)
#mm8.chr4:153178370-155029701 - [2] (4)
#rn3.chr5:171170466-173100353 - [2] (9)
#canFam2.chr5:59136729-61002229 + [2] (95)
#monDom4.chrUn:113503411-113571255 + [0] (1066)
#monDom4.chrUn:113625931-113631959 + [3] (1066)
#monDom4.chr4:380557591-380574536 + [3] (1365)
#monDom4.chr4:380626117-380638649 + [3] (1365)
#monDom4.chr4:380964424-380984297 + [3] (1365)


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

		$line =~ s/(.+)\.(chr.+):(\d+)-(\d+) .+/$1\t$2\t$3\t$4/g;
		# print $line;
        printf($outF "$blockID\t$line\n");

	}
}



close($outF);

print "!!!!!!!!!! DONE !!!!!!!!!!\n";
print "... written in $outputFile\n";


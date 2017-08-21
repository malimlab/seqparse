#!/usr/bin/perl5.16

use strict;

my $file = $ARGV[0];
open (my $fh, "<", $file)
	or die "Cannot open file";

print "Position,Base,A,C,G,T,Coverage,A>A,A>C,A>G,A>T,C>A,C>C,C>G,C>T,G>A,G>C,G>G,G>T,T>A,T>C,T>G,T>T,\n";

while (<$fh>) {
	/HIVss\t(\d+)\t(\w)\t(\d+)\t.*?\tA:(\d+):.*?\tC:(\d+):.*?\tG:(\d+):.*?\tT:(\d+)/;
	my $pos = $1;
	my $base = $2;
	my $readCount = $3;
	my $noA = $4;
	my $noC = $5;
	my $noG = $6;
	my $noT = $7;
	if ($pos > 80) {
		$pos = $pos - 80;
		print "$pos,$base,$noA,$noC,$noG,$noT,$readCount,";
		for ($base) {
			/^C$/ and do {print "0,0,0,0,"; last;};
			/^G$/ and do {print "0,0,0,0,0,0,0,0,"; last;};
			/^T$/ and do {print "0,0,0,0,0,0,0,0,0,0,0,0,"; last;};
		}
		if ($noA > 0) {print sprintf("%.3f", $noA/$readCount);} else {print "0"};
		print ",";
	    if ($noC > 0) {print sprintf("%.3f", $noC/$readCount);} else {print "0"};
		print ",";
		if ($noG > 0) {print sprintf("%.3f", $noG/$readCount);} else {print "0"};
		print ",";
		if ($noT > 0) {print sprintf("%.3f", $noT/$readCount);} else {print "0"};
		print ",";
		for ($base) {
			/^A$/ and do {print "0,0,0,0,0,0,0,0,0,0,0,0"; last;};
			/^C$/ and do {print "0,0,0,0,0,0,0,0,"; last;};
			/^G$/ and do {print "0,0,0,0,"; last;};
		}
		print "\n";
	}
}

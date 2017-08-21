#!/usr/bin/perl5.16

use strict;
use List::Util qw[min max];
use Statistics::Regression;

my @filenames;
my @count;
my @count_adj;
my @count_adj_sum;
my @noalign;
my @readCount_ss;
my $n;
my $max_seqlen;
my $path_base;

foreach my $file (@ARGV) {

	open (my $fh, "<", $file)
	or die "Cannot open file";

	$file =~ /(.*)\/(.*?)_(.*)\.sam/;
	$filenames[$n] = $2;
	$path_base = $1;

	print STDERR "* $2_$3\n";

	while (<$fh>) {
		if (/^.*?\t\d+\tHIVss\t(\d+)\t\d+\t.*?\t.*?\t\d+\t\d+\t(.*?)\t/) {
			my $loc = $1 + length($2) - 81;
			if ($loc > 0) {
				$count[$n][$loc]++;
				$max_seqlen = max($max_seqlen, $loc);
				if ($loc <= 160) {
					$readCount_ss[$n] = $readCount_ss[$n] + 1;
				}
			} else {
				$noalign[$n]++;
			}
		} else {
			$noalign[$n]++;
		}
	}

	close($fh);
	$n++;
}

open (my $fh, ">", "$path_base/../parse_results/readCount.csv");

print $fh ",";
for (my $j = 0; $j < $n; $j++) {
	print $fh "$filenames[$j], "
}
print $fh "\n";

for (my $i = 1; $i <= $max_seqlen; $i++) {
	print $fh "$i, ";
	for (my $j = 0; $j < $n; $j++) {
		if ($count[$j][$i]) {
			print $fh "$count[$j][$i], ";
		} else {
			print $fh "0, ";
		}
	}
	print $fh "\n";
}

print $fh "\nunaligned sequences,";
for (my $j = 0; $j < $n; $j++) {
	print $fh "$noalign[$j], ";
}
print $fh "\n";
close ($fh);

open (my $fh, ">", "$path_base/../parse_results/readCount_normalised.csv");
print $fh ",";
for (my $j = 0; $j < $n; $j++) {
	print $fh "$filenames[$j], "
}
print $fh "\n";

for (my $i = 1; $i <= 160; $i++) {
	print $fh "$i, ";
	for (my $j = 0; $j < $n; $j++) {
		print $fh $count[$j][$i]/$readCount_ss[$j] . ", ";
	}
	print $fh "\n";
}
close ($fh);

print STDERR "Performing length normalisation... ";

my $reg_kwok;
my $reg_nnn;
my @reg_kwok_a;
my $reg_kwok_m;
my $reg_kwok_c;
my @reg_nnn_a;
my $reg_nnn_m;
my $reg_nnn_c;

for (my $j = 0; $j < $n; $j++) {
	if ($filenames[$j] =~ /HTPCon/i) {
		$reg_kwok=Statistics::Regression->new(
			"Kwok", ["Intercept", "Slope"]
		);
		$reg_kwok->include($count[$j][11]/$readCount_ss[$j], [1.0, 11.0]);
		$reg_kwok->include($count[$j][12]/$readCount_ss[$j], [1.0, 12.0]);
		$reg_kwok->include($count[$j][18]/$readCount_ss[$j], [1.0, 18.0]);
		$reg_kwok->include($count[$j][19]/$readCount_ss[$j], [1.0, 19.0]);
		$reg_kwok->include($count[$j][24]/$readCount_ss[$j], [1.0, 24.0]);
		$reg_kwok->include($count[$j][47]/$readCount_ss[$j], [1.0, 47.0]);
		$reg_kwok->include($count[$j][49]/$readCount_ss[$j], [1.0, 49.0]);
		$reg_kwok->include($count[$j][50]/$readCount_ss[$j], [1.0, 50.0]);
		$reg_kwok->include($count[$j][54]/$readCount_ss[$j], [1.0, 54.0]);
		$reg_kwok->include($count[$j][61]/$readCount_ss[$j], [1.0, 61.0]);
		$reg_kwok->include($count[$j][63]/$readCount_ss[$j], [1.0, 63.0]);
		$reg_kwok->include($count[$j][81]/$readCount_ss[$j], [1.0, 81.0]);
		$reg_kwok->include($count[$j][85]/$readCount_ss[$j], [1.0, 85.0]);
		$reg_kwok->include($count[$j][94]/$readCount_ss[$j], [1.0, 94.0]);
		$reg_kwok->include($count[$j][96]/$readCount_ss[$j], [1.0, 96.0]);
		$reg_kwok->include($count[$j][97]/$readCount_ss[$j], [1.0, 97.0]);
		$reg_kwok->include($count[$j][98]/$readCount_ss[$j], [1.0, 98.0]);
		@reg_kwok_a = $reg_kwok->theta();
		$reg_kwok_c = $reg_kwok_a[0];
		$reg_kwok_m = $reg_kwok_a[1];
	}
}

print "\n\nc: $reg_kwok_c\nm: $reg_kwok_m\n";

for (my $j = 0; $j < $n; $j++) {
	if ($reg_kwok) {
		for (my $i = 1; $i <= 160; $i++) {
			$count_adj[$j][$i] = $count[$j][$i] / ($reg_kwok_m * $i + $reg_kwok_c);
			$count_adj_sum[$j] = $count_adj_sum[$j] + $count_adj[$j][$i];
		}
	}
}

open (my $fh, ">", "$path_base/../parse_results/readCount_length_adjusted.csv");
print $fh ",";
for (my $j = 0; $j < $n; $j++) {
	print $fh "$filenames[$j], "
}
print $fh "\n";

for (my $i = 1; $i <= 160; $i++) {
	print $fh "$i,";
	for (my $j = 0; $j < $n; $j++) {
		if ($reg_kwok) {
			print $fh $count_adj[$j][$i]/$count_adj_sum[$j] . ",";
		} else {
			print $fh "NA,";
		}
	}
	print $fh "\n";
}
close ($fh);

print STDERR "done.\n";

if ($reg_kwok) {
	$reg_kwok->print();
	print "\n";
}

if ($reg_nnn) {
	$reg_nnn->print();
	print "\n";
}

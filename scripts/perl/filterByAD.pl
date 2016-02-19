#!/usr/bin/perl
use strict;
use warnings;

(my $pindel_vcf, my $cutoff, my $replacement) = @ARGV;

my $row;
my @genotypes;

open(my $fh, '<', $pindel_vcf) or die "Could not open file '$pindel_vcf' $!";

$pindel_vcf =~ s/\.vcf/.AD_filtered.vcf/; 

open(OUT, ">" , $pindel_vcf);

while($row = <$fh>) {
	if($row =~ m/^#/) {
		print OUT $row;
	} else {
		while($row =~ m/\d\/\d:(\d+),(\d+)/g) {
			
			if($1+$2 < $cutoff) {
				$row =~ s/$&/$replacement/g;
			}
			
		}
		print OUT $row;
	}
}

close(OUT);


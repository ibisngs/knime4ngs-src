#!/usr/bin/perl
#  Copyright (C) 2016 the Knime4NGS contributors.
#  Website: http://ibisngs.github.io/knime4ngs
#  
#  This file is part of the KNIME4NGS KNIME extension.
#  
#  The KNIME4NGS extension is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
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


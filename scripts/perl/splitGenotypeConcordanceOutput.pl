#!/usr/bin/perl -w

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

#title           :splitGenotypeConcordanceOutput.pl
#description     :This script will split the Output of GATK GenotypeConcordance into separate files
#author		 :Maximilian Hastreiter
#date            :20141210
#version         :0.1
#usage		 :perl splitGenotypeConcordanceOutput.pl <INFILE>
#==============================================================================


#Infile
my $INFILE 	= $ARGV[0];
#Outfiles
my @Outfiles 	= ($INFILE."_ProportionsGenotypesComp", $INFILE."_ComparisonCounts", $INFILE."_ProportionsGenotypesEval", $INFILE."_PerSampleSummaryStatistics", $INFILE."_SiteLevelSummaryStatistics");


open(GenotypeConcordance_Outfile, $INFILE);

my $TableCounter = 0;
my $isOpen 	 = 0;

while (<GenotypeConcordance_Outfile>) {
	my $line = $_;

	if($line =~ m/^#/)
	{
		if($isOpen==0){
			if($TableCounter!=0){
				close($curr_Outfile);
				print "Closing $Outfiles[$TableCounter-1]\n"
			}
			print "Opening $Outfiles[$TableCounter]\n";
			open($curr_Outfile, '>', $Outfiles[$TableCounter]);
			$TableCounter++;
			$isOpen = 1;
		}
		print $curr_Outfile $line;
		print $line;
	}else{
		$line =~ s/Non-Reference /Non-Reference_/g;
		$isOpen = 0;
		print $curr_Outfile $line;
		print $line;
	}
		

}


close(GenotypeConcordance_Outfile);

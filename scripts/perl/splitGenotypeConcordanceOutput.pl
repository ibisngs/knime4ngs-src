#!/usr/bin/perl -w
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

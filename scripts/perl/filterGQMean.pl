#!/usr/bin/perl
#title           :filterGQMean.pl
#description     :This script calculates the average GQ value for each site and removes the line if it falls below the specified cutoff. It also removes individual GT entries if GQ value is below GQ Cutoff
#author		 	 :Maximilian Hastreiter
#date            :20150702
#version         :0.2
#usage		 	 :perl filterGQMean.pl <INFILE> <GQ Mean Cutoff> <GQ Cutoff>
#==============================================================================


#Files
my $INFILE 				= $ARGV[0];
my $GQ_MEAN_CUTOFF 		= $ARGV[1];
my $GQ_CUTOFF 			= $ARGV[2];

#Counters
my $processedLines		= 0;
my $filteredLines		= 0;
my $filteredSampleGQ	= 0;


my $missingGTs		= 0;
my $foundGQColumn	= "FALSE";
my $GQIndex			= -1;

open(INVCF,"<", $INFILE) or die "Can not open file $INFILE: $!";

my @INFILE_SPLIT = split(/.vcf/,$INFILE);

open(OUTVCF,">","$INFILE_SPLIT[0].GQFiltered.vcf")  or die "Can not create file $INFILE_SPLIT[0].GQFiltered.vcf: $!";

print "#Running filterGQMean.pl with GQ Cutoff $GQ_CUTOFF and GQ Mean Cutoff $GQ_MEAN_CUTOFF\\n";
print "#Processing $INFILE\\n";
print "#Writing to $INFILE_SPLIT[0].GQFiltered.vcf\\n";


#Iterate File
while(<INVCF>){
	my $line = $_;
	my @fields =split(/\t/,$line);
	if($line =~ m/^#.*/){#Print all Header lines
		print OUTVCF $line;
	}
	
	if(!($line =~ m/^#.*/)){ 	#Variant line	

		chomp($line);
		$processedLines++;

		my @line_fields = split(/\t/,$line);
		my $SumGQ 		= 0;
		my $SampleCount = 1;
		my $outline 	= "";
		
		#Store first 8 cols of orignial line in outline
		$outline 		= join("\t",@line_fields[0..8]); 

		#Search for Index of GQ value
		if($foundGQColumn eq "FALSE"){
			my @SampleInfo = split(/:/,$line_fields[8]);
			for(my $j=-1;$j<scalar(@SampleInfo);$j++){
				
				if($SampleInfo[$j] eq "GQ"){
					$foundGQColumn = "TRUE";
					$GQIndex	   = $j;
					print "#Found GQ value in column $j of FORMAT\\n";
				}
			}
		}


		#Process all samples
		for(my $i=9;$i<scalar(@line_fields);$i++){
			
			if(!($line_fields[$i] =~ m/^\..*/)){						#Skip ./. entries

				my $sample_line = $line_fields[$i];
				my @anno_fields =split(/:/,$sample_line);
				my $GQ = $anno_fields[$GQIndex];
				
				if(!($GQIndex == -1 || $GQ eq ".")){					#Skip entries with . as GQ value or no GQ entry was found
				
						if($GQ>=$GQ_CUTOFF){							#Check if GQ is above the threshold for the specific sample
							
							#Add original sample line to outline
							$outline .= "\t".$sample_line;
							
							#Add GQ value to sum and use it for mean calculation
							$SumGQ = $SumGQ + $GQ;
							$SampleCount++;
							
						}else{#Sample failed the GQ cutoff. GT is set to ./.
							
							$anno_fields[0] 	= "./.";
							my $new_sample_line = join(":",@anno_fields);
							
							$outline .= "\t".$new_sample_line;
							$filteredSampleGQ++;
							#print "############\n";
							#print "$sample_line\n";
							#print "$new_sample_line\n";
							#print "$GQ\n";
							#print "\n";
						}												
				}else{#end skip missing GQ
						$outline .= "\t".$sample_line;	
				}
			}else{#end skip ./.
				$outline .= "\t".$line_fields[$i];
				$missingGTs = $missingGTs + 1 ;	
				#print "$line_fields[$i]\n";		
	
			}
			
		}#end for
			
		my $MeanGQ = $SumGQ/$SampleCount;
		if($MeanGQ >= $GQ_MEAN_CUTOFF){
			print OUTVCF $outline."\n";
		}else{

			$filteredLines++;
		
		}

	    }#end if variant line

}#end while
					     
print "#Processed $processedLines variants\\n";	
print "#Found $missingGTs missing genotypes\\n";
print "#Found $filteredSampleGQ genotypes with GQ below $GQ_CUTOFF\\n";	
print "#Filtered $filteredLines variants due to GQMean < $GQ_MEAN_CUTOFF\\n";			  



close(INVCF);
close(OUTVCF);

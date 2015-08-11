package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

/**
 * 
 * @author Charu
 * 1) parse the thunder file
 * 2) get samples and save the information
 * 3) write comments to vcf
 * 4) check contigs available in the thunder output and write to vcf
 * 5) write the references according to the contigs
 * 6) write the main information to vcf
 */
public class ThunderToVCF {

	private Set<String> contig;
	private ArrayList<String> samples = new ArrayList<String>();
	private ArrayList<String> thunder = new ArrayList<String>();


	
//	private void setContig (String s) {
//		this.contig = s;
//	}
//	
//	private void setReference (String s) {
//		this.reference = s;
//	}
	private void getThunderInfo (String thunderFile) throws IOException{
		ArrayList<String> tmp_contig = new ArrayList<String>();
		
		File input = new File (thunderFile);
		try{
			Scanner sc = new Scanner(input, "UTF-8");
			String[] header = sc.nextLine().trim().split("\t");
			if ((header.length - 6) %3 != 0) { // chr pos ref al1/al2 posterior freq1
				System.err.println("Thunder format wrong");
				System.exit(0);
			}
			for (int i = 6 ; i < header.length; i = i+3) {
//				System.out.println(header[i]);
				this.samples.add(header[i].substring(header[i].lastIndexOf("/")).split("_")[1]);
			}
//			for (String s : this.samples)
//				System.out.println(s);
			
	        while (sc.hasNextLine()) {
	        	String currentLine = sc.nextLine().trim();
	        	if(currentLine == null )
		            throw new IOException( thunderFile + " File is empty " ); 
	        	String[] tmp = currentLine.split("\t");
//	        	System.out.println(tmp.length);
	        	String[] sampleInfo = new String[(tmp.length - 6)/3];
	        	int j = 0;
	        	for (int i = 6; i < tmp.length; i = i+3) {
	        		if (tmp[i].equals("0"))
	        			sampleInfo[j] = "0/0:";
	        		else if (tmp[i+1].equals("0"))
	        			sampleInfo[j] = "0/1:";
	        		else if (tmp[i+2].equals("0"))
	        			sampleInfo[j] = "1/1:";
	        		else {
	        			int a = Integer.parseInt(tmp[i]),
	        					b = Integer.parseInt(tmp[i+1]),
	        					c = Integer.parseInt(tmp[i+2]);
	        			if (a < b && a < c)
	        				sampleInfo[j] = "0/0";
	        			else if (b < a && b < c)
	        				sampleInfo[j] = "0/1";
	        			else if (c < a && c < b)
	        				sampleInfo[j] = "1/1";
	        		}
	        		j++;
	        	}
	        	
	        	tmp_contig.add(tmp[0]); // add the chromosome into the list while reading line
	    		NumberFormat formatter = new DecimalFormat("#0.000");
	        	String info = tmp[0] + "\t" //chromosome
	        				+ tmp[1] + "\t" //position
	        				+ "." + "\t" // ID
	        				+ tmp[2] + "\t" // ref
	        				+ tmp[3].substring(2) + "\t" // alt
	        				+ tmp[4] + "\t" // QUAL/posterior
	        				+ "." + "\t" // Filter
	        				+ "AF=" + formatter.format(1.0 - Double.parseDouble(tmp[5])) + "\t" // INFO
	        				+ "GT:PL"; // FORMAT
	        	int k = 0;
	        	for (int i = 6; i < tmp.length; i = i+3) {
	        		info += "\t" + sampleInfo[k] + tmp[i] + "," + tmp[i+1] + "," + tmp[i+2];
	        		k++;
	        	}
	        	info += "\n";
	        	
//	        	System.out.print(info);
	        	this.thunder.add(info);
	        }
		}
		catch (FileNotFoundException e) {
	        e.printStackTrace();
	    }
		this.contig = new HashSet<String>(tmp_contig); //fetching the unique choromosomes
//		for(String s : this.contig)
//			System.out.println(s);
	}

	private void writeComments (BufferedWriter bw) {
		FastaCheck fc = new FastaCheck();
		
		try { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
			String s = "##fileformat=VCFv4.1" + "\n" +
						"##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">" + "\n" +
						"##FILTER=<ID=LowQual,Description=\"Low quality\">" + "\n" +
						"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" + "\n" +
						"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">" + "\n" +
						"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" + "\n" +
						"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" + "\n" +
						"##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP observed within the GVCF block\">" + "\n" +
						"##FORMAT=<ID=PGT,Number=1,Type=String,Description=\"Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another\">" + "\n" +
						"##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group\">" + "\n" +
						"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">" + "\n" +
						"##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.\">" + "\n" +
						"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">" + "\n" +
						"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">" + "\n" +
						"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">" + "\n" +
						"##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">" + "\n" +
						"##INFO=<ID=CCC,Number=1,Type=Integer,Description=\"Number of called chromosomes\">" + "\n" +
						"##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases\">" + "\n" +
						"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">" + "\n" +
						"##INFO=<ID=DS,Number=0,Type=Flag,Description=\"Were any of the samples downsampled?\">" + "\n" +
						"##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">" + "\n" +
						"##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">" + "\n" +
						"##INFO=<ID=GQ_MEAN,Number=1,Type=Float,Description=\"Mean of all GQ values\">" + "\n" +
						"##INFO=<ID=GQ_STDDEV,Number=1,Type=Float,Description=\"Standard deviation of all GQ values\">" + "\n" +
						"##INFO=<ID=HWP,Number=1,Type=Float,Description=\"P value from test of Hardy Weinberg Equilibrium\">" + "\n" +
						"##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description=\"Consistency of the site with at most two segregating haplotypes\">" + "\n" +
						"##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">" + "\n" +
						"##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed\">" + "\n" +
						"##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed\">" + "\n" +
						"##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">" + "\n" +
						"##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total Mapping Quality Zero Reads\">" + "\n" +
						"##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">" + "\n" +
						"##INFO=<ID=NCC,Number=1,Type=Integer,Description=\"Number of no-called samples\">" + "\n" +
						"##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">" + "\n" +
						"##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">" + "\n" +
						"##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">" + "\n" +
						"##Thunder_call_command=very complicated" + "\n";
			bw.write(s);
			
			for (String chr: this.contig) {
				int length = 0;
				for (int i = 0; i < fc.chr_length.length; i++) {
					if (chr.equals(fc.chr_length[i].getId()))
						length = fc.chr_length[i].getLength();
				}
					
				bw.write("##contig=<ID=" + chr + ",length=" + length + ">"+ "\n"); //##contig=<ID=chr21,length=48129895>
				bw.write("##reference=file:///home/ibis/tanzeem.haque/Documents/hg19/" + chr + ".fa" + "\n"); //##reference=file:///home/ibis/tanzeem.haque/Documents/hg19/chr21.fa
			}
			
			bw.write("#CHROM" + "\t" + "POS" + "\t" + "ID" + "\t" + "REF" + "\t" + "ALT" 
			+ "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" + "FORMAT");
			for (String indiv : this.samples)
				bw.write("\t" + indiv);
			bw.write("\n");
			
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}
	}
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		ThunderToVCF t2vcf = new ThunderToVCF();
		t2vcf.getThunderInfo("fin21_75RL.chr21"); // get the infos parsed, get the contigs
		BufferedWriter bw = new BufferedWriter(new FileWriter("fin21_75RL.chr21.vcf"));
		t2vcf.writeComments(bw);
		for (String s : t2vcf.thunder)
			bw.write(s);
		bw.close();

	}

}

package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;

/**
 * @author tanzeem.haque
 * This class checks how many and which chromosomes are included in the (multi)fasta file
 * and their respective lengths are also saved hier
 *
 */
public class FastaCheck {
	
	public final static class HG_19 {
		
		private String id;
		private int length
		;
		
		HG_19(String id, int length) {
			setId(id);
			setLength(length);
		}

		public String getId() {
			return id;
		}

		public void setId(String id) {
			this.id = id;
		}

		public int getLength() {
			return length;
		}

		public void setLength(int length) {
			this.length = length;
		}
		
		public String toString() {
			return getId() + "\t" + (getLength()/10000000) + 1 + "\n";
		}

	}

	/**
	 * chr_length stores the length of the respective chromosome
	 * samtools faidx Homo_sapiens_assembly18.fasta 
	 */
	public final HG_19[] chr_length = {new HG_19("chr1",249250621), new HG_19("chr2",243199373), new HG_19("chr3",198022430), new HG_19("chr4",191154276),
			 new HG_19("chr5",180915260), new HG_19("chr6",171115067), new HG_19("chr7",159138663), new HG_19("chr8",146364022), 
			 new HG_19("chr9",141213431), new HG_19("chr10",135534747), new HG_19("chr11",135006516), new HG_19("chr12",133851895), 
			 new HG_19("chr13",115169878), new HG_19("chr14",107349540), new HG_19("chr15",102531392), new HG_19("chr16",90354753), 
			 new HG_19("chr17",81195210), new HG_19("chr18",78077248), new HG_19("chr19",59128983), new HG_19("chr20",63025520), 
			 new HG_19("chr21",/*50*/48129895), new HG_19("chr22",/*50*/51304566), new HG_19("chrX",155270560), new HG_19("chrY",59373566)};
//	
//	final HG_19[] chr_length = {new HG_19("seq1", 447)};// {new HG_19("seq1",10263), new HG_19("seq2",3056)};
	
	/**
	 * input_chr_length stores the chr_length thingy according to the chromosomes in the (multi)fasta input
	 * 
	 * Also for ID_Chunksize for recordfile ids_chunks.txt
	 */
	ArrayList<String> input_chr_length = new ArrayList<String>(25); // can be maximum 24 chromosomes
	/**
	 * so that we can know from the very beginning how many mutations we need
	 */
	long input_length = 0;
	/**
	 * choose the F0 or F1 and M0 or M1 in the very beginninng
	 */
	private String[] parentalChromatids;
	
	public FastaCheck() {
		
	}

	FastaCheck(String file) throws IOException {
		calculate_length(file);
		choose_parentalChromatids();
	}
	
	
	/**
	 * @throws IOException 
	 * find out the length of the whole input
	 * @throws  
	 */
	void calculate_length(String file) throws IOException {
		File input = new File (file);
		try {
			Scanner sc = new Scanner(input, "UTF-8");
			
	        while (sc.hasNextLine()) {
	        	String currentLine = sc.nextLine().trim();
	        	if(currentLine == null )
		            throw new IOException( file + " File is empty " );     
//		        if(currentLine.charAt(0) != '>')
//		            throw new IOException("Header of the file " + file + " must be '>'");

	        	if(currentLine.length() > 0 && currentLine.charAt(0) == '>'){
	        		currentLine = currentLine.substring(1).trim();//.replaceAll("\\s+", "_");
//	        		System.out.println("ID: " + currentLine);
	        		int tmp_l = getLength(currentLine);
	        		this.input_chr_length.add(currentLine + "\t" + tmp_l);
	        		this.input_length += (long) tmp_l;	 
//	        		System.out.println(tmp_l);
//	        		FrostRunner.memory();

	        	} 
//	        	else 
//	        		continue;
	        		
	        }
	        
	    } 
	    catch (FileNotFoundException e) {
	        e.printStackTrace();
	    }
	}

	private int getLength(String currentLine) {
		// TODO Auto-generated method stub
		int l = 0;
		for (HG_19 hg : chr_length) {
			if (currentLine.equals(hg.getId()))
				l = hg.getLength();
		}		
		return l;

	}
	
	/**
	 * @param fileID : the 4 files from parents
	 * @param id
	 * @return the two files to get the mutated sequences to inherit
	 */
	String[] choose_parentalChromatids() {
		// TODO Auto-generated method stub
		String[] readFile_for_child = new String[2];
		int calc_prob = probability();

		switch(calc_prob) {

		case 0: //F0M0
			setParentalChromatids("M0", "F0");
			readFile_for_child[0]="_M_0.fa";
			readFile_for_child[1]="_F_0.fa";
			break;
		case 1:  //F0M1
			setParentalChromatids("M0", "F1");
			readFile_for_child[0]="_M_0.fa";
			readFile_for_child[1]="_F_1.fa";
			break;
		case 2: //F1M0
			setParentalChromatids("M1", "F0");
			readFile_for_child[0]="_M_1.fa";
			readFile_for_child[1]="_F_0.fa";
			break;
		case 3: //F1M1
			setParentalChromatids("M1", "F1");
			readFile_for_child[0]="_M_1.fa";
			readFile_for_child[1]="_F_1.fa";
			break;

		}
//		System.out.println("I CAME IN TO CHOOSE FILE");
//		for (String s : readFile_for_child) 
//			System.out.println(s);

		return readFile_for_child;
	}
	
	/**
	 * Calculates the probability for a random combination of the alleles
	 * random combination of alleles = 2^23 (for human) = 8388608
	 * probability of a diploid combination: 2^23!/((2^23-2)!*2!) = 3.5184367894528*10^13,  100/... = 2.84*10^(-12)
	 */
	private int probability () {
		Random rand = new Random();
		int i = rand.nextInt(4); // 4 combinations of chromosomes
		return i;

	}
	
	//a1, a2, id
		protected String[] getParentalChromatids() {
			return parentalChromatids;
		}

		/**
		 * @param al1
		 * @param al2
		 * Id basically not necessary since the combination is already decided from the very beginning
		 */
		protected void setParentalChromatids (String al1, String al2){ //, String id) {
			this.parentalChromatids = new String[2];
			this.parentalChromatids[0] = al1;
			this.parentalChromatids[1] = al2;

		}
}

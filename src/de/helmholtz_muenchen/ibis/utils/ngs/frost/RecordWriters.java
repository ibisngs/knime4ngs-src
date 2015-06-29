package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;



/**
 * @author tanzeem.haque
 * Class to write all info outputs
 *
 */

/**
 * @author tanzeem.haque
 * Class to write all info outputs
 *
 */

public class RecordWriters {

	private FastaReader fs;
	private InputScanner in;
	private TrioSimulator trio;
	private ArrayList<VCF_info> vcf = new ArrayList<>();
	private ArrayList<Integer> duplicates = new ArrayList<>();
	
	public RecordWriters(FastaReader fs, InputScanner in, TrioSimulator trio) {
		// TODO Auto-generated constructor stub
		setFs(fs);
		setIn(in);
		setTrio(trio);
	}
	
	public RecordWriters() {
		
	}
	
	protected void write_simple_string(BufferedWriter bw, String s) {
		// TODO Auto-generated method stub
		
		try { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
			bw.write(s);
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}
//		File file = new File(fileName);
//		try {
//			BufferedWriter bw = new BufferedWriter( new FileWriter(file));
//			bw.write(s + "\n");	
//			bw.close();
//		} catch (Exception e) {
//			System.err.println("Error: " + e.getMessage());
//			e.printStackTrace();
//		}		
	}
	
	protected void write_InputData(BufferedWriter bw, InputData iData/*, int n*/) {
		// TODO Auto-generated method stub
		
		try { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
			for (int i = 0; i < iData.getPositions()./*get(n).*/size(); i++) {
				bw.write(iData.getId() + "\t" + iData.getPositions()./*get(n).*/get(i) + "\n");
				
			}
			
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}		
	}
	
	protected void write_strand(BufferedWriter bw_strand, String vcf_File) {
		// TODO Auto-generated method stub
		ArrayList<Integer> mutations = new ArrayList<>();
		File vcf = new File (vcf_File);
		try {
			Scanner sc = new Scanner(vcf, "UTF-8");
			while (sc.hasNextLine()) {
				String s = sc.nextLine();
				if (s.startsWith("#"))
					continue;
				String[] col = s.split("\t");
				int pos = Integer.parseInt(col[1]);
				mutations.add(pos);
			}
		}
		catch (FileNotFoundException e) {
	        e.printStackTrace();
	    }
		ArrayList<String> strandMap = prepareStrandFile(mutations);
		
		String strand = "";
		for(int i = 0; i < strandMap.size(); i++) {
			strand += strandMap.get(i) + "\n"; 
		}
		write_simple_string(bw_strand, strand);
		
	}
	private ArrayList<String> prepareStrandFile(ArrayList<Integer> mutations){
		// TODO Auto-generated method stub
		ArrayList<String> strandMap = new ArrayList<>(mutations.size());
		
		for (int n : mutations) {
			for(int i = 0; i < in.getStrandMap().size(); i++) {
				if (n >= in.getStrandMap().get(i).getStart() && n <= in.getStrandMap().get(i).getEnd()){
					String strand = (n) + "\t" + in.getStrandMap().get(i).getStream();
					strandMap.add(strand);
					break;
				}			
			}
		}
		return strandMap;
	}
	
	protected void write_vcf(BufferedWriter bw_vcf, String parentFile, String childFile) {
		File child_tmp = new File (childFile), parents_tmp = new File (parentFile);
		String m = "", f = ""; // father allele(F0/F1), mother allele (M0/M1)
		try {
			Scanner sc = new Scanner(child_tmp, "UTF-8");
			
	        while (sc.hasNextLine()) {
	        	String s = sc.nextLine();
//	        	System.out.println(s);
	        	String [] col = s.trim().split("\t");
	        	int pos = Integer.parseInt(col[1]);
	        	m = col[6];
	        	f = col[7];
	        	if (col[3].equals("*")) {
//	        		System.out.println("*************************");
	        		col[3] = col[3].replace("*", "");
	        		pos = pos - 1;
	        		
	        	}
	        	this.vcf.add(new VCF_info(pos /*position*/, 
	        			col[0]/*.replaceAll("^[0-9]_", "")/*id*/
	        			+"\t"+pos
	        			+"\t"+"." /*ID as in vcf format*/	
	        			+"\t"+col[2]/*ref allele*/
	        			+"\t"+col[3] /*alt allele*/
	        			+"\t"+"0.0" /*QUAL as in vcf format*/
	        			+"\t"+"." /* FILTER as in vcf format */
	        			+"\t"+"AN=6" /* INFO as in vcf format*/
	        			+"\t"+"GT" /* FORMAT as in vcf format */
	        			+"\t"+"0|0" /*mother genotype*/
	        			+"\t"+"0|0" /*father genotype*/
	        			+"\t"+col[4].replace("/", "|")/*child genotype*/));
	        }
			truncateFiles(child_tmp);

	        sc = new Scanner(parents_tmp, "UTF-8");
	        while(sc.hasNextLine()) {
	        	String [] col = sc.nextLine().trim().split("\t");
	        	String c0 = "", c1 ="";
	        	c0 = (m.equals("M0"))?col[4]/*mother genotype*/.split("/")[0]:col[4].split("/")[1];
	        	c1 = (f.equals("F0"))?col[5]/*father genotype*/.split("/")[0]:col[5].split("/")[1];
	        	int pos = Integer.parseInt(col[1]);
	        	if (col[3].equals("*")) {
//	        		System.out.println("*************************");
	        		col[3] = col[3].replace("*", "");
	        		pos = pos - 1;
	        		
	        	}
	        	this.vcf.add(new VCF_info(pos, 
	        			col[0]/*.replaceAll("^[0-9]_", "")/*id*/
	    	        			+"\t"+col[1]
	    	    	        	+"\t"+"." /*ID as in vcf format*/	
	    	        			+"\t"+col[2]/*ref allele*/
	    	        			+"\t"+col[3]/*alt allele*/
	    	        			+"\t"+"0.0" /*QUAL as in vcf format*/
	    	    	        	+"\t"+"." /* FILTER as in vcf format */
	    	    	        	+"\t"+"AN=6" /* INFO as in vcf format*/
	    	    	        	+"\t"+"GT" /* FORMAT as in vcf format */
	    	        			+"\t"+col[4].replace("/", "|")/*mother genotype*/
	    	    	        	+"\t"+col[5].replace("/", "|")/*father genotype*/
	    	    	        	+"\t"+c0+"|"+c1/*child genotype*/));
	        }
			truncateFiles(parents_tmp);

	        Collections.sort(this.vcf);
	        
	        removeDuplicates(this.duplicates);
	        
		}
		catch (FileNotFoundException e) {
	        e.printStackTrace();
	    }
		try  { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
//			System.out.println("It should write the vcf!");
			
			for (VCF_info v : this.vcf) {
//				System.out.println(v.getContent());
				bw_vcf.write(v.getContent()+"\n");
			}
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}
		
	}

	private void truncateFiles(File tmpFile) {
	// TODO Auto-generated method stub
		File txtFile = new File (tmpFile.toString().replace(".tmp", ".txt"));
		ArrayList<String> data = new ArrayList<>();
		
		/**
		 * save the info from the recordfiles (.tmp) in another file
		 */
		try {
			Scanner sc = new Scanner(tmpFile, "UTF-8");		
	        while (sc.hasNextLine()) {
	        	String s = sc.nextLine();
	        	data.add(s);
	        }
		}
	    catch (FileNotFoundException e) {
		       e.printStackTrace();
		}
		/**
		 * truncate .tmp
		 */
	    try {
			FileChannel outChan = new FileOutputStream(tmpFile, true).getChannel();
			outChan.truncate(0);
		    outChan.close();

		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		/**
		 * write the saved data from the .tmps to .txt
		 */
		try /*(BufferedWriter pw = new BufferedWriter(new FileWriter(txtFile, true), 10000000)) {*/(PrintWriter pw = new PrintWriter(new FileOutputStream(txtFile, true)))
		{
			for (String s : data) {
				pw.write(s + "\n");
				
			}			
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}	
		
	}

	protected void unphase (BufferedWriter bw, String vcfFile) {
		ArrayList<String> data = new ArrayList<>();
		File phased = new File (vcfFile);
//		File unphased = new File (vcfFile.replace("phased", "unphased"));
		try {
			Scanner sc = new Scanner(phased, "UTF-8");
			
	        while (sc.hasNextLine()) {
	        	String s = sc.nextLine();
	        	s = s.replaceAll("\\|", "/");
	    		s = s.replaceAll("1/0", "0/1");
	    		s = s.replaceAll("2/1", "1/2");
	    		s = s.replaceAll("2/0", "0/2");
	    		data.add(s);
//	        	System.out.println(s);
	        	
	        }
	          
		}
		catch (FileNotFoundException e) {
	        e.printStackTrace();
	    }
		try { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
//			System.out.println("It should write the vcf!");
    		for (String s : data) 
    			bw.write(s+"\n");
		
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}
		
	}
	
	private void removeDuplicates(ArrayList<Integer> duplicates) {
		// TODO Auto-generated method stub

		String s = (duplicates.size() > 0)? "Found double mutation in child =>"
				:"No double mutation in child";		
		FrostRunner.createLog(FrostRunner.bw_log , s);
//		for (int i : duplicates)
//			System.out.println(i + "");

		for (int dup_position : duplicates) {
			for (int i = 0; i < this.vcf.size(); i++) {
				if (dup_position == this.vcf.get(i).getPosition()) {
					/**
					 * the first entry is child entry due to the data structure
					 */
					/**
					 * saving the parental mutation to add to the denovo mut //record from parents file
					 */
					int c0 = (FrostRunner.parental_chromatids[0].equals("M0"))?
							Integer.parseInt(this.vcf.get(i+1).getContent().split("\t")[9]/*mother genotype*/.split("\\|")[0])
							:Integer.parseInt(this.vcf.get(i+1).getContent().split("\t")[9].split("\\|")[1]);
		        	int c1 = (FrostRunner.parental_chromatids[1].equals("F0"))?
		        			Integer.parseInt(this.vcf.get(i+1).getContent().split("\t")[10]/*father genotype*/.split("\\|")[0])
							:Integer.parseInt(this.vcf.get(i+1).getContent().split("\t")[10].split("\\|")[1]);
		        	/**
		             * denovo situation
		             */
		            String [] c = this.vcf.get(i).getContent().split("\t")[11].split("\\|"); //record from child file
		            String denovo = this.vcf.get(i).getContent().split("\t")[4], parent_mut = this.vcf.get(i+1).getContent().split("\t")[4];

		            /**
		             * KEEP this system out block to check bug
		             */
		            FrostRunner.createLog(FrostRunner.bw_log,"Child from parent: " + c0 + "|" + c1);
		            FrostRunner.createLog(FrostRunner.bw_log,"Denovo situation: " + c[0] + "|" + c[1]);
		            FrostRunner.createLog(FrostRunner.bw_log,"Child allele: " + denovo + ";" + "\t" + "Parent allele: " + parent_mut);

		            if(!(denovo.equals(parent_mut))
		            		|| (denovo.equals(parent_mut) && (c0 == 0 || c1 == 0))) { 
		            /**
		             * child nucleotide is not the same as parents, that means an authentic denovo
		             */
//		            	System.out.println("I ADDED");
		            	c0 += Integer.parseInt(c[0]);//m
		            	c1 += Integer.parseInt(c[1]);//f
		            }
		            /**
		             * col[0].replaceAll("^[0-9]_", "")//id
	    	        			+"\t"+col[1]
	    	    	        	+"\t"+"." //ID as in vcf format	
	    	        			+"\t"+col[2]//ref allele
	    	        			+"\t"+col[3]//alt allele
	    	        			+"\t"+"0.0" //QUAL as in vcf format
	    	    	        	+"\t"+"." // FILTER as in vcf format 
	    	    	        	+"\t"+"AN=6" // INFO as in vcf format
	    	    	        	+"\t"+"GT" // FORMAT as in vcf format 
	    	        			+"\t"+col[4].replace("/", "|")//mother genotype
	    	    	        	+"\t"+col[5].replace("/", "|")//father genotype
	    	    	        	+"\t"+c0+"|"+c1//child genotype));
		             */
		            String new_content = this.vcf.get(i).getContent().split("\t")[0]+"\t"
    						+ this.vcf.get(i).getContent().split("\t")[1]+"\t"
    						+ this.vcf.get(i).getContent().split("\t")[2]+"\t"
    						+ this.vcf.get(i+1).getContent().split("\t")[3]+"\t"
    						+ this.vcf.get(i).getContent().split("\t")[4]+"\t"/*mutated base of the child -> twice*/
    						+ this.vcf.get(i).getContent().split("\t")[5]+"\t"
    	    				+ this.vcf.get(i).getContent().split("\t")[6]+"\t"
    	    				+ this.vcf.get(i).getContent().split("\t")[7]+"\t"
    	    				+ this.vcf.get(i).getContent().split("\t")[8]+"\t"
    						/**
    						 * BEWARE: STANDARDIZED M F C INSTEAD OF F M C 
    						 *(Already changed in the class Mutation (getParentString))
    						 * C always M/F
    						 */
    						+ this.vcf.get(i+1).getContent().split("\t")[9]+"\t"/*mother genotype originally from parents file*/
    						+ this.vcf.get(i+1).getContent().split("\t")[10]+"\t"/*father genotype originally from parents file*/
    						+ c0 + "|" + c1;
		            FrostRunner.createLog(FrostRunner.bw_log, new_content);
		            this.vcf.set(i, new VCF_info(this.vcf.get(i).getPosition(), new_content));
//    				System.out.println("Before removing: " + tmp.get(i).getContent());
		            this.vcf.remove(i+1);
//    				System.out.println("After removing: " + tmp.get(i).getContent());
		            break;
				}
			}
		}
	}

	/**
	 * @return the fs
	 */
	public FastaReader getFs() {
		return fs;
	}

	/**
	 * @param fs the fs to set
	 */
	public void setFs(FastaReader fs) {
		this.fs = fs;
	}

	/**
	 * @return the in
	 */
	public InputScanner getIn() {
		return in;
	}

	/**
	 * @param in the in to set
	 */
	public void setIn(InputScanner in) {
		this.in = in;
	}

	/**
	 * @return the family
	 */
	public TrioSimulator getTrio() {
		return trio;
	}

	/**
	 * @param trio the family to set
	 */
	public void setTrio(TrioSimulator trio) {
		this.trio = trio;
	}
	
	/**
	 * Nested Class to create a cpmparable object for vcf
	 */

	private class VCF_info implements Comparable<VCF_info>{
		private int position;
		private String content;
		
		private VCF_info(int position, String content) {
			setPosition(position);
			setContent(content);
		}

		public int getPosition() {
			return position;
		}

		public void setPosition(int position) {
			this.position = position;
		}

		public String getContent() {
			return content;
		}

		public void setContent(String content) {
			this.content = content;
		}

		@Override
		public int compareTo(VCF_info o) {
			// TODO Auto-generated method stub
			int comp = this.getPosition()-o.getPosition();
			if (comp==0 && this.getContent().split("\t")[0].equals(o.getContent().split("\t")[0])){
				duplicates.add(o.getPosition());
//				System.out.println("SAME");

			}

			return comp;
		}
	}

	
}
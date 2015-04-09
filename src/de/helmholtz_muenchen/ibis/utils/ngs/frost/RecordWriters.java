package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
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
	private Trio_Simulator trio;
	private ArrayList<VCF_info> vcf = new ArrayList<>();
	private ArrayList<Integer> duplicates = new ArrayList<>();
	
	public RecordWriters(FastaReader fs, InputScanner in, Trio_Simulator trio) {
		// TODO Auto-generated constructor stub
		setFs(fs);
		setIn(in);
		setTrio(trio);
	}
	
	public RecordWriters() {
		
	}
	
	protected void write_simple_string(String fileName, String s) {
		// TODO Auto-generated method stub
		
		try (PrintWriter pw = new PrintWriter(new FileOutputStream(fileName, true))) { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
			pw.write(s);
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
	
	protected void write_InputData(String fileName, InputData iData, int n) {
		// TODO Auto-generated method stub
		
		try (PrintWriter pw = new PrintWriter(new FileOutputStream(fileName, true))) { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
			for (int i = 0; i < iData.getPositions().get(n).size(); i++) {
				pw.write(iData.getId() + "\t" + iData.getPositions().get(n).get(i) + "\n");
				
			}
			
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}		
	}
	
//	protected void write_vcf(String fileName, ArrayList<String> ID, 
//			ArrayList<Integer> POS, ArrayList<String> prettyMutation) {
//		
//		// TODO Auto-generated method stub
//		try (PrintWriter pw = new PrintWriter(new FileOutputStream(fileName, true))) { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
//			for (int i = 0; i < prettyMutation.size(); i++) {
//				pw.write(ID.get(i) + "\t" + POS.get(i) + "\t" + prettyMutation.get(i) + "\n");
//			}
//		} catch (Exception e) {
//			System.err.println("Error: " + e.getMessage());
//			e.printStackTrace();
//		}		
//	}
	
	protected void write_vcf(String outFile, String parentFile, String childFile) {
		File child = new File (childFile), parents = new File (parentFile);
		String f = "", m = ""; // father allele(F0/F1), mother allele (M0/M1)
		try {
			Scanner sc = new Scanner(child, "UTF-8");
			
	        while (sc.hasNextLine()) {
	        	String s = sc.nextLine();
//	        	System.out.println(s);
	        	String [] col = s.trim().split("\t");
	        	f = col[6];
	        	m = col[7];
	        	this.vcf.add(new VCF_info(Integer.parseInt(col[1]) /*position*/, 
	        			col[0]/*id*/
	        			+"\t"+col[1]
	        			+"\t"+col[2]/*ref allele*/
	        			+"\t"+col[3] /*alt allele*/
	        			+"\t"+"0/0" /*father genotype*/
	        			+"\t"+"0/0" /*mother genotype*/
	        			+"\t"+col[4]/*child genotype*/));
	        }
	        
	        sc = new Scanner(parents, "UTF-8");
	        while(sc.hasNextLine()) {
	        	String [] col = sc.nextLine().trim().split("\t");
	        	String c0 = "", c1 ="";
	        	c0 = (m.equals("M0"))?col[5]/*mother genotype*/.split("/")[0]:col[5].split("/")[1];
	        	c1 = (f.equals("F0"))?col[4]/*father genotype*/.split("/")[0]:col[4].split("/")[1];
	        	this.vcf.add(new VCF_info(Integer.parseInt(col[1]), 
	        			col[0]/*id*/
	    	        			+"\t"+col[1]
	    	        			+"\t"+col[2]/*ref allele*/
	    	        			+"\t"+col[3]/*alt allele*/
	    	        			+"\t"+col[4]/*father genotype*/
	    	    	        	+"\t"+col[5]/*mother genotype*/
	    	    	        	+"\t"+c0+"/"+c1/*child genotype*/));
	        }
	        Collections.sort(this.vcf);
	        removeDuplicates(this.duplicates);
	        
		}
		catch (FileNotFoundException e) {
	        e.printStackTrace();
	    }
		try (PrintWriter pw = new PrintWriter(new FileOutputStream(outFile))) { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
//			System.out.println("It should write the vcf!");
			for (VCF_info v : this.vcf)
				pw.write(v.getContent()+"\n");
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}
	}

	private void removeDuplicates(ArrayList<Integer> duplicates) {
		// TODO Auto-generated method stub
		ArrayList<VCF_info> tmp = new ArrayList<>(this.vcf.size());
		tmp = this.vcf;
		int duplicate_idx = 0;
		for (int i : duplicates) {
			i = i - duplicate_idx;
			/**
			 * saving the parental mutation to add to the denovo mut
			 */
			int c0 = (FrostRunner.parental_chromatids[1].equals("M0"))?
					Integer.parseInt(tmp.get(i).getContent().split("\t")[5]/*mother genotype*/.split("/")[0])
					:Integer.parseInt(tmp.get(i).getContent().split("\t")[5].split("/")[1]);
        	int c1 = (FrostRunner.parental_chromatids[0].equals("F0"))?
        			Integer.parseInt(tmp.get(i).getContent().split("\t")[4]/*father genotype*/.split("/")[0])
					:Integer.parseInt(tmp.get(i).getContent().split("\t")[4].split("/")[1]);
        	
        	String [] c = tmp.get(i-1).getContent().split("\t")[6].split("/");//f/m
        	c0 += Integer.parseInt(c[1]);//m
        	c1 += Integer.parseInt(c[0]);//f
        	String new_content = tmp.get(i).getContent().split("\t")[0]+"\t"
        						+ tmp.get(i).getContent().split("\t")[1]+"\t"
        						+ tmp.get(i).getContent().split("\t")[2]+"\t"
        						+ tmp.get(i-1).getContent().split("\t")[3]+"\t"/*mutated base of the child -> twice*/
        						+ tmp.get(i).getContent().split("\t")[4]+"\t"
        						+ tmp.get(i).getContent().split("\t")[5]+"\t"
        						+ c0 + "/" + c1;
        	tmp.set(i, new VCF_info(tmp.get(i).getPosition(), new_content));
//        	System.out.println("Before removing: " + tmp.get(i).getContent());
        	tmp.remove(i-1);
//        	System.out.println("After removing: " + tmp.get(i-1).getContent());
//        	System.out.println("After removing: " + tmp.get(i).getContent());
        	duplicate_idx++;
		}
		this.vcf= tmp;
		
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
	public Trio_Simulator getTrio() {
		return trio;
	}

	/**
	 * @param trio the family to set
	 */
	public void setTrio(Trio_Simulator trio) {
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
			if (comp==0){
				duplicates.add(vcf.indexOf(o));
//				System.out.println("SAME");
//				System.out.println(this.getContent());//child
//				System.out.println(o.getContent());//parents
//				
//	        			
//	        	int i = vcf.indexOf(o);
//	        	int j = vcf.indexOf(this);
//	        	System.out.println(i + "\t" + j);
			}

			return comp;
		}
	}
}
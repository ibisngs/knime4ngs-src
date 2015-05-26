package de.helmholtz_muenchen.ibis.utils.ngs.frost;

//package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

public class VCF_DataCollection {
	
	private static String indiv1;
	private static String indiv2;
	private static String indiv3;
		
	protected static ArrayList<TrioData> extractVcf (String anyVcfFile) throws IOException {
		ArrayList<TrioData> tdal = new ArrayList<>();
		File input = new File (anyVcfFile);
		
		try {
			Scanner sc = new Scanner(input, "UTF-8");
			
	        while (sc.hasNextLine()) {
        		boolean fileOk = false;
	        	String currentLine = sc.nextLine().trim();
	        	if(currentLine == null )
		            throw new IOException( anyVcfFile + " File is empty " );     

	        	/**
	        	 * skip the vcf info part
	        	 */
	        	if(currentLine.length() > 0 && currentLine.startsWith("##"))
	        		continue;	    
	        	
	        	else {
	        		/**
	        		 * Check the clumn length at #CHROM
	        		 */
	        		if(currentLine.length() > 0 && currentLine.startsWith("#CHROM")){
		        		String[] cols = currentLine.split("\t");
		        		if (cols.length == 12) {
		        			VCF_DataCollection.indiv1 = cols[9];//usually M, can be otherwise
		        			VCF_DataCollection.indiv2 = cols[10];//usually F
		        			VCF_DataCollection.indiv3 = cols[11];//Usually C
		        			fileOk = true;
		        			continue;
		        		}
	        		}
	        		if(fileOk) {
	        			boolean gtOk = false;
	        			String[] cols = currentLine.split("\t");
	        			String pos = cols[1], ref = cols[3], alt = cols[4],
	        					mom = trioOrder(cols)[0], pop = trioOrder(cols)[1], kid = trioOrder(cols)[2];
	        			if (mom.length()==3 && pop.length()==3 && kid.length()==3)
		        			gtOk = true;
	        			
	        			if (gtOk) {
	        				int i = Integer.parseInt(pos);
	        				Genotype f = new Genotype(Integer.parseInt(pop.substring(0, 1)), Integer.parseInt(pop.substring(2, 3)));
	        				Genotype m = new Genotype(Integer.parseInt(mom.substring(0, 1)), Integer.parseInt(mom.substring(2, 3)));
	        				Genotype c = new Genotype(Integer.parseInt(kid.substring(0, 1)), Integer.parseInt(kid.substring(2, 3)));

	        				TrioData td = new TrioData(i, /*position*/
	        						new Allelic_Change(ref, alt),/*Ref and alt*/
	        						new TrioGenotype(f,m,c) /*F, M, C*/);
	        				
	        				tdal.add(td);
	        			}
	        			else {
	        				System.err.println("Something wrong with column and GT; GT length not 3");
	        				System.exit(0);
	        			}
	        		}
	        	}
	        }
	    } 
	    catch (FileNotFoundException e) {
	        e.printStackTrace();
	    }
		
		return tdal;
	}

	private static String[] trioOrder(String[] cols) {
		String mom = "", pop = "", kid = "";
		
		switch(VCF_DataCollection.indiv1) {
		case "M":
			mom = cols[9];
		case "F":
			pop = cols[9];
		case "C":
			kid = cols[9];
			break;
		}
		switch(VCF_DataCollection.indiv2) {
		case "M":
			mom = cols[10];
		case "F":
			pop = cols[10];
		case "C":
			kid = cols[10];
			break;
		}
		switch(VCF_DataCollection.indiv3) {
		case "M":
			mom = cols[11];
		case "F":
			pop = cols[11];
		case "C":
			kid = cols[11];
			break;
		}
		String[] order = {mom, pop, kid};
		return order;
	} 

	protected static boolean checkContig (String anyVcfFile1, String anyVcfFile2) throws IOException {
		
		boolean contigSame = false;
		String[] contig1 = VCF_DataCollection.contig(anyVcfFile1);
		String[] contig2 = VCF_DataCollection.contig(anyVcfFile2);
		
		if (contig1[0].equals(contig2[0])/*ID*/ && contig1[1].equals(contig2[1])/*length*/)
			contigSame = true;
			
		return contigSame;
	}
	
	private static String[] contig (String anyVcfFile) throws IOException {
		
		File input = new File (anyVcfFile);
		String contig = "", ID = "", length = "";
		try {
			Scanner sc = new Scanner(input, "UTF-8");
			
	        while (sc.hasNextLine()) {
	        	String currentLine = sc.nextLine().trim();
	        	if(currentLine == null )
		            throw new IOException( anyVcfFile + " File is empty " );     

	        	if(currentLine.length() > 0 && currentLine.startsWith("##contig")){
	        		contig = currentLine; //"##contig=<ID=chr21,length=48129895>";
	        		ID = contig.substring(contig.indexOf("<")+"ID=".length()+1, contig.indexOf(","));
	        		length = contig.substring(contig.indexOf(",")+"length=".length()+1, contig.indexOf(">"));
	        		break;
	        	} 
	        }        
	    } 
	    catch (FileNotFoundException e) {
	        e.printStackTrace();
	    }
		String[] con = {ID, length};
		return con;
	}
	
}

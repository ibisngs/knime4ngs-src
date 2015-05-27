package de.helmholtz_muenchen.ibis.ngs.filtergatkvariants;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import org.knime.core.node.InvalidSettingsException;

public class FilterGATKVariants {
	
	/**
	 * Removes all variants of the Pindel VCF that already exist in GATK VCF File
	 * @param GATK_Indel_File
	 * @param Pindel_File
	 * @throws InvalidSettingsException 
	 */
	public static String filter(String GATK_Indel_File, String Pindel_File) throws InvalidSettingsException{
		
		String OUT_FILE = Pindel_File.replace(".vcf", ".filteredGATKVariants.vcf");
		HashMap<String,String> GATK_VCF = new HashMap<String, String>();
		String GATK_VCF_HEADER			= "";
		String PINDEL_VCF_HEADER		= "";
		//Load GATK Indels
		GATK_VCF_HEADER = loadGATKVariants(GATK_VCF, GATK_Indel_File,GATK_VCF_HEADER);
		
		//Go through Pindel VCF and filter Pindel<->GATK Duplicates
		try(BufferedReader br = new BufferedReader(new FileReader(Pindel_File))) {
	        String line;
			
			PrintWriter writer = new PrintWriter(OUT_FILE, "UTF-8");
	        
	        while ((line = br.readLine()) != null) {
	        	
	        	line = line.replace(System.getProperty("line.separator"), "");
	        	
	        	if(line.startsWith("#")){
	        		writer.println(line);
	        		//Load Pindel Header
	        		if(line.startsWith("#CHROM")){
	        			String[] fields = line.split("\t");
	        			for(int i=9;i<fields.length;i++){
	        				if(i==9){
	        					PINDEL_VCF_HEADER = fields[i].split("_")[0];
	        				}else{
	        					PINDEL_VCF_HEADER += "\t"+fields[i].split("_")[0];
	        				}
	        			}
	        			
		        		//Compare Header Lines
		        		if(!PINDEL_VCF_HEADER.equals(GATK_VCF_HEADER)){
		        			throw new InvalidSettingsException("Sample Names do not match!\n"+PINDEL_VCF_HEADER+"\n vs. \n"+GATK_VCF_HEADER);
		        		}
	        		}
	        		

	        		
	        	}else{
		        	String[] fields = line.split("\t");
		        	String key 		= fields[0]+"\t"+fields[1];
		        	if(GATK_VCF.containsKey(key)){
		        		
		        		//Get both GT Strings
		        		String[] GATK_GTs  = GATK_VCF.get(key).split("\t");
		        		String[] PindelGTs = getGTString(fields[3], fields[4], line).split("\t");
		        		
		        		//Outstring
		        		String Filtered_GTs = fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+fields[4]+"\t"+fields[5]+"\t"+fields[6]+"\t"+fields[7]+"\t"+fields[8];
		        		
		        		//Compare and remove dups
		        		for(int i=0;i<GATK_GTs.length;i++){
		        			
		        			if(GATK_GTs[i].equals(PindelGTs[i])){		//Pindel Variant already in GATK, therefore we remove it
		        				Filtered_GTs += "\t"+"./.";			
		        			}else{
		        				Filtered_GTs += "\t"+fields[9+i];
		        			}
		        			
		        		}
		        		
		        		writer.println(Filtered_GTs);
		        		
		        	}else{
		        		writer.println(line);
		        	}
	        	}	        	
	        }
	        
	        br.close();
	        writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return OUT_FILE;
	}
	
	
	/**
	 * Load Variants into HashMap
	 * @param GATK_VARIANTS
	 * @param GATK_Indel_File
	 */
	private static String loadGATKVariants(HashMap<String, String> GATK_VARIANTS, String GATK_Indel_File, String GATK_VCF_HEADER){

		try(BufferedReader br = new BufferedReader(new FileReader(GATK_Indel_File))) {
	        String line;

	        while ((line = br.readLine()) != null) {
	        	if(!(line.startsWith("#"))){
	        		line = line.replace(System.getProperty("line.separator"), "");
	        		String[] fields = line.split("\t");
	   
	        		String GTString = getGTString(fields[3], fields[4], line);     		
	        		GATK_VARIANTS.put(fields[0]+"\t"+fields[1], GTString);

	        	}else{
	        		//Store Header Line of GATK VCF
	        		if(line.startsWith("#CHROM")){
	        			String[] fields = line.split("\t");
	        			GATK_VCF_HEADER = "";
	        			for(int i=9;i<fields.length;i++){
	        				if(i==9){
	        					GATK_VCF_HEADER = fields[i];
	        				}else{
	        					GATK_VCF_HEADER += "\t"+fields[i];
	        				}
	        			}
	        		}
	        		
	        	}
	        }
	        
	        br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return GATK_VCF_HEADER;
	}
	
	
	/**
	 * Convert A	T	0/1	1/1... to A/T	T/T....
	 * @param REF
	 * @param ALT
	 * @param line
	 * @return
	 */
	private static String getGTString(String REF, String ALT, String line){
		
		String GTString 	= "";
		String[] fields 	= line.split("\t");
		String AlleleString = REF+","+ALT;
		String[] Alleles 	= AlleleString.split(",");
		boolean first 		= true;
		
		for(int i=9;i<fields.length;i++){
		
			String currGT = fields[i].split(":")[0];
				
			if(!currGT.equals("./.")){
				
				int GT_1 = Integer.parseInt(currGT.split("/")[0]);
				int GT_2 = Integer.parseInt(currGT.split("/")[1]);
				
				String Allele_1 = Alleles[GT_1];
				String Allele_2 = Alleles[GT_2];
			
				if(first){
					GTString = Allele_1+"/"+Allele_2;
					first	 = false;
				}else{
					GTString += "\t"+Allele_1+"/"+Allele_2;
				}
			
			}else{
				
				if(first){
					GTString = "./.";
					first	 = false;
				}else{
					GTString += "\t"+"./.";
				}
			}
		}

		return GTString;
	}
	
	
}

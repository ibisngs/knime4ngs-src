package de.helmholtz_muenchen.ibis.ngs.fillmissinggenotypes;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import org.knime.core.node.InvalidSettingsException;

/**
 * @author hastreiter
 * Fills missing GTs when merging single-calling Pindel results. Variants with DP<20 will be left as ./. 
 */


public class FillPindelGenotypes {
	
/**
 * 
 * @param VCF_INFILE
 * @param PATH_DATA_DIR
 * @param SUFFIX
 * @throws IOException
 * @throws InvalidSettingsException 
 */
	protected static String fillGTs(String VCF_INFILE, String PATH_DATA_DIR, String SUFFIX) throws IOException, InvalidSettingsException{
		
		/**
		 * Find Coverage Files
		 */
		LinkedList<String> CoverageFiles = new LinkedList<String>();
		de.helmholtz_muenchen.ibis.utils.ngs.FileSearch.search(PATH_DATA_DIR, SUFFIX, CoverageFiles);
		
		HashMap<String,BufferedReader> FileReaders = new HashMap<String,BufferedReader>();
		HashMap<String,String> LastLine = new HashMap<String,String>();
		
		
		//Get BufferedReader for each of the Coverage Files
		for(String File: CoverageFiles){
			
	    	String[] FilePathParts  = File.split("/");
	    	String FileName 		= FilePathParts[FilePathParts.length-1];
	    	String SampleName 		= FileName.split("_rev")[0];
	    	System.out.println(SampleName+"\t"+File);
	    	FileReaders.put(SampleName, new BufferedReader(new FileReader(File)));
		}
		
		
		//Open the VCF File
		String OUT_FILE = VCF_INFILE.replace(".vcf", ".FilledGTs.vcf");
		PrintWriter writer = new PrintWriter(OUT_FILE, "UTF-8");
		String[] Header = null;
		
		//Go through vcf file
	    try(BufferedReader br = new BufferedReader(new FileReader(VCF_INFILE))) {
	        
	    	String line;
	        while ((line = br.readLine()) != null) {

	        	line = line.replace(System.getProperty("line.separator"), "");
	        	
	        	if(line.startsWith("#")){ //Header Line
	        		if(line.startsWith("#CHROM")){
	        			Header = line.split("\t");
	        		
		        		//Check if all Coverage files are available
		        		for(int i=9;i<Header.length;i++){
		        			if(!FileReaders.containsKey(Header[i])){
		        				throw new InvalidSettingsException("Sample ID "+Header[i]+" does not match any of the detected coverage files! Something seems to be wrong with the sample ids");
		        			}
		        		}
	        		}
	        		//Write header to new file#
	        		writer.println(line);
	        	
	        	}else{//Variant line
	        		
	        		//Split the line and get relevant information
	        		String[] fields = line.split("\t");
	        		String chr_variant = fields[0];
	        		int pos_number_variant = Integer.parseInt(fields[1]);
	        		
	        		
	        		/**
	        		 * Create Line with filled GTs
	        		 */
	        		String outline="";
	        		/**
	        		 * Iterate for each sample
	        		 */
	        		for(int i = 0;i<fields.length;i++){
	        			
	        			if(i<9){
	        				//Keep Info fields
	        				if(i==0){
	        					outline=fields[i];
	        				}else{
	        					outline+="\t"+fields[i];
	        				}
	        				
	        			}else{//Check GT for each sample
	        				
		        			String currSampleID = Header[i];
		        			String currEntry 	= fields[i];
		        			
		        			if(currEntry.equals("./.") || currEntry.equals("0/0")){		//Only check if Genotype equals ./. or 0/0
		  	
			        			String Sample_Line;	//Holds the coverage value
			        			while((Sample_Line = FileReaders.get(currSampleID).readLine()) != null){
			        				
			        				String lastLine = ""; //Last processed line
			        				//Load the previously processed line
			        				if(LastLine.containsKey(currSampleID)){
			        					lastLine = LastLine.get(currSampleID);
			        				}
			        				
			        				if(!Sample_Line.startsWith("Locus")){//Skip header line

			        					String[] Sample_LineFields  = Sample_Line.split("\t");
				        				String[] ChrPos 			= Sample_LineFields[0].split(":");
				        				String Chr_curr 			= ChrPos[0];
				        				
				        				//Current Position
				        				int Pos_number_curr = Integer.parseInt(ChrPos[1]);

				        				//If Position of VCF File was found in Coverage file
				        				if(Sample_Line.startsWith(chr_variant+":"+pos_number_variant)){
				        					
				        					LastLine.put(currSampleID, Sample_Line);	//Store the line
				        					
				        					//Correct Position found
				        					int Coverage = Integer.parseInt(Sample_LineFields[1]);
				        					if(Coverage>=20){
				        						outline+="\t"+"0/0:"+Coverage+",0";
				        					}else{
				        						//Keep old entry
				        						outline+="\t"+"./.";
				        					}
				        					
				        					//Next sample
				        					break;
				        					

				        				}else if((!chr_variant.equals(Chr_curr)) || pos_number_variant < Pos_number_curr){
				        					//Already passed the region...No entry found
				        				
				        					//Check last line since there might be multiple VCF entries with the same position after merging!
				        					if(lastLine.startsWith(chr_variant+":"+pos_number_variant)){
				        						//Correct Position found
					        					int Coverage = Integer.parseInt(Sample_LineFields[1]);
					        					if(Coverage>=20){
					        						outline+="\t"+"0/0:"+Coverage+",0";
					        					}else{
					        						//Keep old entry
					        						outline+="\t"+"./.";
					        					}
				        					}else{
				        						LastLine.put(currSampleID, Sample_Line);	//Store the line
				        						//Keep old entry
				        						outline+="\t"+"./.";
				        					}
				        					//Next sample
				        					break;
				        			
				        				}else{
				        					
				        					//Keep on searching
				        					
				        				}
			        				}			    	
			        			}
			        			
		        			}else{
		        				//Keep old entry
		        				outline+="\t"+currEntry;
		        			}
	        			}
	        		}
	        		//Write new vcf line
	        		writer.println(outline);
	        	}
	        	
	        }
	        
	        //Close all FHs
		    Iterator<String> iterator_FR = FileReaders.keySet().iterator();  
		    while (iterator_FR.hasNext()) {
		    		String key = iterator_FR.next();
		    		FileReaders.get(key).close();
		    }
	        
	        writer.close();
	        br.close();
	    }
		
		return OUT_FILE;
		
	}
	
	
	
}

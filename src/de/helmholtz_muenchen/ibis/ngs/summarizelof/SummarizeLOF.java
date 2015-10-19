package de.helmholtz_muenchen.ibis.ngs.summarizelof;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import de.helmholtz_muenchen.ibis.utils.ngs.PEDFile;

public class SummarizeLOF {

	
	public SummarizeLOF(){
		
	}
	
	/**
	 * 
	 * @param PED_FILE
	 * @param VCF_FILE
	 * @return
	 * @throws UnsupportedEncodingException 
	 * @throws FileNotFoundException 
	 */
	@SuppressWarnings("unchecked")
	public String[] getLOFs(String PED_FILE, String VCF_FILE) throws FileNotFoundException, UnsupportedEncodingException{
		
		//Read PED File
		HashMap<String,String> PED = new PEDFile(PED_FILE).getPED();
		
		//Find LOFs in VCF
		HashMap<String,HashMap<String,String>> LOF_GENES = findLOFs(VCF_FILE);
		
		//Get Summary Stats
		HashMap<String,ArrayList<String>>[] LOFSummaries = getLOFSummary(LOF_GENES);
		HashMap<String,ArrayList<String>> SAMPLE_SUMMARY= LOFSummaries[0];
		HashMap<String,ArrayList<String>> GENE_SUMMARY= LOFSummaries[1];
		
		
		//2 more outfiles for LOF genes
		String CHILD_LOFS		 					= VCF_FILE.replace("vcf", "childLOFs.csv");
	    String SampleHOM_ParentsHET_FULLLOF_FILE	= VCF_FILE.replace("vcf", "SampleHOM_ParentsHET_FULLLOF.csv");
	    PrintWriter writer_CHILD_LOFS = new PrintWriter(CHILD_LOFS, "UTF-8");
	    PrintWriter writer_SampleHOM_ParentsHET_FULLLOF = new PrintWriter(SampleHOM_ParentsHET_FULLLOF_FILE, "UTF-8");
		
		//Get Hom-Het Stats
		Iterator<String> iterator = PED.keySet().iterator();  
	       
	    while (iterator.hasNext()) { 
		       String SampleID 	= iterator.next().toString(); 
		       String Parents = PED.get(SampleID);
		       String PaternalID = Parents.split("\t")[0];
		       String MaternalID = Parents.split("\t")[1];

		       if(LOF_GENES.containsKey(SampleID)){
		    	
		    	   HashMap<String,String> Sample_LOF_GENES	 = LOF_GENES.get(SampleID);
		    	   
		    	   if(LOF_GENES.containsKey(PaternalID) && LOF_GENES.containsKey(MaternalID)){
		    		   HashMap<String,String> Paternal_LOF_GENES = LOF_GENES.get(PaternalID);
				       HashMap<String,String> Maternal_LOF_GENES = LOF_GENES.get(MaternalID);
				          
				       Iterator<String> inner_iterator = Sample_LOF_GENES.keySet().iterator(); 
				       int SampleHOM_ParentsALL				= 0;
				       int SampleHOM_ParentsHET				= 0;
				       int SampleHOM_ParentsALL_FULLLOF		= 0;
				       int SampleHOM_ParentsHET_FULLLOF		= 0;
				       
				       while (inner_iterator.hasNext()) {
				    	
					       String Gene 	= inner_iterator.next().toString(); 
					       String Info = Sample_LOF_GENES.get(Gene);
					       String SampleHomHet = Info.split(";")[0];
					       String SampleFullLOF = Info.split(";")[3];
				    	   
					       writer_CHILD_LOFS.println(Gene+"\t"+SampleID);
					       
					       if(Paternal_LOF_GENES.containsKey(Gene) && Maternal_LOF_GENES.containsKey(Gene) && SampleHomHet.equals("hom")){
					    	
					    	   SampleHOM_ParentsALL++;	
					    	   
					    	   String PaternalHomHet = Paternal_LOF_GENES.get(Gene).split(";")[0];
					    	   String MaternalHomHet = Maternal_LOF_GENES.get(Gene).split(";")[0];
					    	   String PaternalFullLOF = Paternal_LOF_GENES.get(Gene).split(";")[3];
					    	   String MaternalFullLOF = Maternal_LOF_GENES.get(Gene).split(";")[3];
					    	   
					    	   if(SampleFullLOF.equals("TRUE") && PaternalFullLOF.equals("TRUE") && MaternalFullLOF.equals("TRUE")){
					    		   SampleHOM_ParentsALL_FULLLOF++;
					    	   }
					    	   
					    	   
					    	   if(PaternalHomHet.equals("het") && MaternalHomHet.equals("het")){
					    		   SampleHOM_ParentsHET++;
					    		   
						    	   if(SampleFullLOF.equals("TRUE") && PaternalFullLOF.equals("TRUE") && MaternalFullLOF.equals("TRUE")){
						    		   SampleHOM_ParentsHET_FULLLOF++;
						    		   writer_SampleHOM_ParentsHET_FULLLOF.println(Gene+"\t"+SampleID);
						    		   
						    		   System.out.println(Gene+"\t"+SampleID);
						    	   }
					    		   
					    	   } 
					       }
				       }
				       
				       //Get stored values
				       ArrayList<String> SAMPLE_VALUE = SAMPLE_SUMMARY.get(SampleID);
				       ArrayList<String> PATERNAL_VALUE = SAMPLE_SUMMARY.get(PaternalID);
				       ArrayList<String> MATERNAL_VALUE = SAMPLE_SUMMARY.get(MaternalID);
				       SAMPLE_VALUE.add(""+SampleHOM_ParentsALL);
				       SAMPLE_VALUE.add(""+SampleHOM_ParentsALL_FULLLOF);
				       SAMPLE_VALUE.add(""+SampleHOM_ParentsHET);
				       SAMPLE_VALUE.add(""+SampleHOM_ParentsHET_FULLLOF);
				       PATERNAL_VALUE.add("-");
				       PATERNAL_VALUE.add("-");
				       PATERNAL_VALUE.add("-");
				       PATERNAL_VALUE.add("-");
				       MATERNAL_VALUE.add("-");
				       MATERNAL_VALUE.add("-");
				       
				       SAMPLE_SUMMARY.put(SampleID,SAMPLE_VALUE);
				       SAMPLE_SUMMARY.put(PaternalID,PATERNAL_VALUE);
				       SAMPLE_SUMMARY.put(MaternalID,MATERNAL_VALUE);
		    	   }
		       }
	    }
		
	    /**
	     * Write Outfiles
	     */
	    String SAMPLE_SUMMARY_OUT 			= VCF_FILE.replace("vcf", "lof_SampleSummary.csv");
	    String SAMPLE_SUMMARY_HEADER 		= "SAMPLE\tstartOverlap\tendOverlap\tprematureStop\tremovedStop\tspliceOverlap\tinsertionFS\tdeletionFS\tNumberOfLOFs\tNumberOfFullLOFs\tSampleHOM_ParentsALL\tSampleHOM_ParentsALL_FullLOF\tSampleHOM_ParentsHET\tSampleHOM_ParentsHETFullLOF";
	    String GENE_SUMMARY_OUT 			= VCF_FILE.replace("vcf", "lof_GeneSummary.csv");
	    String GENE_SUMMARY_HEADER 			= "Gene\tNumberOfAffectedSamples\tDeleteMe";
		printMap(SAMPLE_SUMMARY,SAMPLE_SUMMARY_OUT,SAMPLE_SUMMARY_HEADER);
		printMap(GENE_SUMMARY,GENE_SUMMARY_OUT,GENE_SUMMARY_HEADER);	
		
		//Close File writers
	    writer_CHILD_LOFS.close();
	    writer_SampleHOM_ParentsHET_FULLLOF.close();
		
		return new String[]{SAMPLE_SUMMARY_OUT,GENE_SUMMARY_OUT};
	}
	

	
	/**
	 * Reads a VCF File and stores all LOF Mutations in a HashMap SampleID-->{LOFGene1;[het,hom];GT;Type;FullLOF	LOFGene2;...}
	 * @param VCF_FILE
	 * @return
	 */
	private HashMap<String, HashMap<String,String>> findLOFs(String VCF_FILE ){
		
		HashMap<String,HashMap<String,String>> LOF_GENES = new HashMap<String,HashMap<String,String>>();
		
		try{
	  		  FileInputStream fstream = new FileInputStream(VCF_FILE);
	  		  // Get the object of DataInputStream
	  		  DataInputStream in = new DataInputStream(fstream);
	  		  BufferedReader br = new BufferedReader(new InputStreamReader(in));
	  		  String strLine;

	  		  String[] HEADER = {};
	  		  
	  		  while ((strLine = br.readLine()) != null)   {
	  			  
	  			  if(strLine.startsWith("#")){//Commentline
	  				  if(strLine.startsWith("#CHR")){//Headerline
	  					  HEADER = strLine.split("\t");
	  				  }
	  			  }else{//Variant
	  				  /**Get Annotation**/
	  				  String[] fields 		= strLine.split("\t");
	  				  String INFO_FIELD 	= fields[7];
	  				  String[] ANNOTATIONS 	= INFO_FIELD.split(";");
	  				  String VAT_ANNOTATION = "NA";
	  				  String LOF_ALLELES = "";
	  				  String LOF_GENE = "";
	  				  String LOF_TYPE = "";
	  				  String LOF_TRANSCRIPT = "";
	  				  
	  				  /**Get VAT Annotation**/
	  				  for(String a : ANNOTATIONS){
	  					  if(a.startsWith("VA=")){
	  							VAT_ANNOTATION = a;  
	  				
	  							break;
	  					  }
	  				  }
	  				   				  
	  				  /**If VAT Annotation is available**/
	  				  if(!VAT_ANNOTATION.equals("NA")){	  					  
/**
* VA=1:MUC20:ENSG00000176945.11:+:nonsynonymous:1/6:MUC20-011:ENST00000423938.1:382_89_30_P->Q
1:MUC20:ENSG00000176945.11:+:synonymous:5/6:MUC20-002:ENST00000436408.1:2169_1855_619_R->R:MUC20-201:ENST00000320736.6:1614_1342_448_R->R:MUC20-001:ENST00000447234.2:2127_1855_619_R->R:MUC20-004:ENST00000445522.2:2022_1750_584_R->R:MUC20-202:ENST00000381954.5:1560_1288_430_R->R
2:MUC20:ENSG00000176945.11:+:nonsynonymous:1/6:MUC20-011:ENST00000423938.1:382_89_30_P->L
2:MUC20:ENSG00000176945.11:+:prematureStop:5/6:MUC20-002:ENST00000436408.1:2169_1855_619_R->*:MUC20-201:ENST00000320736.6:1614_1342_448_R->*:MUC20-001:ENST00000447234.2:2127_1855_619_R->*:MUC20-004:ENST00000445522.2:2022_1750_584_R->*:MUC20-202:ENST00000381954.5:1560_1288_430_R->*
*/
	  					  
	  					  //Precheck: Check if Annotation line contains LOF Keywords
	  					  if(VAT_ANNOTATION.contains("prematureStop") || VAT_ANNOTATION.contains("removedStop") || VAT_ANNOTATION.contains("spliceOverlap") || VAT_ANNOTATION.contains("insertionFS") || VAT_ANNOTATION.contains("deletionFS") || VAT_ANNOTATION.contains("startOverlap") || VAT_ANNOTATION.contains("endOverlap")){
	  						 //Remove 'VA=' Tag
	  						 VAT_ANNOTATION = VAT_ANNOTATION.replaceFirst("VA=", ""); 
	  						 
		  					  String[] VAT_SUB_ANNOTATION 	= VAT_ANNOTATION.split(",");
		  					  
		  					  
		  					  for(String annotation : VAT_SUB_ANNOTATION){
		  						  
		  						 
		  						  
		  						  String[] VAT_FIELDS = annotation.split(":");
			  					  String AlleleNumber 	= VAT_FIELDS[0];
			  					  String Gene 			= VAT_FIELDS[1];
			  					  String Type 			= VAT_FIELDS[4];
			  					  
			  					  if(Type.contains("|")){ 
			  						  Type = Type.split("|")[0];
			  					  }		  					
			  					  
			  					  String Transcript		= VAT_FIELDS[5];
			  					   
			  					  if(Type.equals("prematureStop") || Type.equals("removedStop") || Type.equals("spliceOverlap") || Type.equals("insertionFS") || Type.equals("deletionFS") || Type.equals("startOverlap") || Type.equals("endOverlap")){
			  						  LOF_ALLELES+=AlleleNumber+";"; 
			  						  LOF_GENE+=Gene+";"; 
			  						  LOF_TYPE+=Type+";";
			  						  LOF_TRANSCRIPT+=Transcript+";";
			  						  System.out.println(LOF_TYPE);
			  					  }
		  					  }
		  					  String [] LOF_ALLELES_FIELDS = LOF_ALLELES.split(";");
		  					  String [] LOF_TYPE_FIELDS = LOF_TYPE.split(";");
		  					  String [] LOF_TRANSCRIPT_FIELDS = LOF_TRANSCRIPT.split(";");
		  					  
		  					  //Now search for affected individuals
		  					  for(int i = 9; i < fields.length;i++){
		  						 
		  						  String SampleID 			= HEADER[i]; 
		  						  String IndividualEntry 	= fields[i];
		  						  String GT 				= IndividualEntry.split(":")[0];
		  						  
		  						  //Workaround for Pindel GTs "."
		  						  if(GT.equals(".")){
		  							  GT = "./.";
		  						  }
		  						  String[] GT_fields 		= GT.split("[/|]");
		  						  String Type 				= "";
		  						  String FullLOF			= "";
		  							  
		  						  //Index of the GT Alleles in the LOF_TYPE_FIELDS Array; -1 indicates missing value
		  						  int IndexFirstAllele  = contains(GT_fields[0], LOF_ALLELES_FIELDS);
		  						  int IndexSecondAllele = contains(GT_fields[1], LOF_ALLELES_FIELDS);
		  						  		  						  
		  						  if(IndexFirstAllele!=-1 && IndexSecondAllele!=-1){	//If Sample contains two LOF Alleles
		  							 
		  							  //Get the correct Mutation Type
		  							  if(LOF_TYPE_FIELDS[IndexFirstAllele].equals(LOF_TYPE_FIELDS[IndexSecondAllele])){
		  								  Type = LOF_TYPE_FIELDS[IndexFirstAllele];
		  							  }else{
		  								  Type = LOF_TYPE_FIELDS[IndexFirstAllele]+"|"+LOF_TYPE_FIELDS[IndexSecondAllele];
		  							  }
		  							  
		  							  //Check Transcripts
		  							  String[] AffectedTranscriptsFirstAllele = LOF_TRANSCRIPT_FIELDS[IndexFirstAllele].split("/");
		  							  String[] AffectedTranscriptsSecondAllele = LOF_TRANSCRIPT_FIELDS[IndexSecondAllele].split("/");
		  							  if(AffectedTranscriptsFirstAllele[0].equals(AffectedTranscriptsFirstAllele[1]) && AffectedTranscriptsSecondAllele[0].equals(AffectedTranscriptsSecondAllele[1])){
		  								FullLOF = "TRUE";
		  							  }else{
		  								FullLOF = "FALSE";  
		  							  }
		  							  
		  							  
		  							  if(LOF_GENES.containsKey(SampleID)){
		  								  
		  								  String Gene = LOF_GENE.split(";")[0];				//Get Gene of first Transcript entry
		  								  LOF_GENES.get(SampleID).put(Gene, "hom;"+GT+";"+Type+";"+FullLOF);
		  								 
		  							  }else{
		  								  
		  								  String Gene = LOF_GENE.split(";")[0];				//Get Gene of first Transcript entry
		  								  HashMap<String,String> tmp = new HashMap<String,String>();
		  								  tmp.put(Gene, "hom;"+GT+";"+Type+";"+FullLOF);
		  								  LOF_GENES.put(SampleID, tmp);
		  							  }
		  							  
		  							  
		  						  }else if (IndexFirstAllele!=-1 || IndexSecondAllele!=-1){ //If Sample contains at least one LOF Allele
		  							  
		  							//Get the correct Mutation Type and Transcripts
		  							  if(IndexFirstAllele!=-1){
		  								Type = LOF_TYPE_FIELDS[IndexFirstAllele]; 
		  								String[] AffectedTranscriptsFirstAllele = LOF_TRANSCRIPT_FIELDS[IndexFirstAllele].split("/");
			  							  if(AffectedTranscriptsFirstAllele[0].equals(AffectedTranscriptsFirstAllele[1])){
				  								FullLOF = "TRUE";
				  							  }else{
				  								FullLOF = "FALSE";  
				  							  }
		  							  }else{
		  								Type = LOF_TYPE_FIELDS[IndexSecondAllele];  
			  							String[] AffectedTranscriptsSecondAllele = LOF_TRANSCRIPT_FIELDS[IndexSecondAllele].split("/");
			  							  if(AffectedTranscriptsSecondAllele[0].equals(AffectedTranscriptsSecondAllele[1])){
				  								FullLOF = "TRUE";
				  							  }else{
				  								FullLOF = "FALSE";  
				  							  }
		  							  }
		  							  
		  							  /**
		  							   * Add to HashMap
		  							   */
		  							  if(LOF_GENES.containsKey(SampleID)){
		  								  
		  								  String Gene = LOF_GENE.split(";")[0];				//Get Gene of first Transcript entry
		  								  LOF_GENES.get(SampleID).put(Gene, "het;"+GT+";"+Type+";"+FullLOF);
		  								 
		  							  }else{
		  								  
		  								  String Gene = LOF_GENE.split(";")[0];				//Get Gene of first Transcript entry
		  								  HashMap<String,String> tmp = new HashMap<String,String>();
		  								  tmp.put(Gene, "het;"+GT+";"+Type+";"+FullLOF);
		  								  LOF_GENES.put(SampleID, tmp);
		  							  }
		  							  
		  						  }
		  						  
		  					  }
		  					  
	  					  }	  					  
	  				  }
	  			  }
	  		  }
	  		br.close();
		}catch (Exception e){//Catch exception if any
			System.err.println("Error: " + e.getMessage());
		}	
		
		return LOF_GENES;
		
	}
	
	/**
	 * Produces two HashMaps(Gene and Sample LOF Summaries)
	 * @param LOF_GENES
	 * @return
	 */
	@SuppressWarnings("rawtypes")
	private HashMap[] getLOFSummary(HashMap<String,HashMap<String,String>> LOF_GENES){
		
		HashMap<String, ArrayList<String>> SAMPLE_SUMMARY = new HashMap<String,ArrayList<String>>();
		HashMap<String, ArrayList<String>> GENE_SUMMARY = new HashMap<String,ArrayList<String>>();
		
		Iterator<String> iterator = LOF_GENES.keySet().iterator();  
	       
	    while (iterator.hasNext()) { 
		       String SampleID 	= iterator.next(); 
		       HashMap<String,String> SampleLOFs = LOF_GENES.get(SampleID);
		       
		       int pmStop 	= 0;
		       int rmStop 	= 0;
		       int spliceO 	= 0;
		       
		       int insFS	= 0;
		       int delFS	= 0;
		       
		       int stOv		= 0;
		       int endOv	= 0;
		       
		       int fullLOFs = 0;
		       
		       int NumberOfLOFs = SampleLOFs.size();
		       
		       Iterator<String> inner_iterator = SampleLOFs.keySet().iterator(); 
		       while (inner_iterator.hasNext()) {
		    	   
		    	   String Gene = inner_iterator.next();
		    	   String Info = SampleLOFs.get(Gene);
		    	   String Type = Info.split(";")[2];
		    	   String FullLOF = Info.split(";")[3];
		    	   
		    	   //Gene Count
		    	   if(GENE_SUMMARY.containsKey(Gene)){
		    		   
		    		   ArrayList<String> OldValue = GENE_SUMMARY.get(Gene);
		    		   String LastSampleIDInserted = OldValue.get(1);

		    		   if(!SampleID.equals(LastSampleIDInserted)){
			    		   int GeneCount = Integer.parseInt(OldValue.get(0));
			    		   GeneCount++;
			    		   GENE_SUMMARY.put(Gene, new ArrayList<String>(Arrays.asList(""+GeneCount,SampleID)));
		    		   }
		    		   
		    	   }else{
		    		   GENE_SUMMARY.put(Gene, new ArrayList<String>(Arrays.asList("1",SampleID)));
		    	   }
		    	   
		    	   //Mutation Counts
		    	   if(Type.equals("prematureStop")){
		    		   pmStop++;
		    	   }else if(Type.equals("removedStop")){
		    		   rmStop++;
		    	   }else if(Type.equals("spliceOverlap")){
		    		   spliceO++;
		    	   }else if(Type.equals("insertionFS")){
		    		   insFS++;
		    	   }else if(Type.equals("deletionFS")){
		    		   delFS++;
		    	   }else if(Type.equals("startOverlap")){
		    		   stOv++;
		    	   }else if(Type.equals("endOverlap")){
		    		   endOv++;
		    	   }else{
		    		   System.out.println(Type); 
		    	   }
		    	   
		    	   if(FullLOF.equals("TRUE")){
		    		   fullLOFs++;
		    	   }
		    	   
		       }
		       //Sample Summary with Counts
		       SAMPLE_SUMMARY.put(SampleID, new ArrayList<String>(Arrays.asList(""+stOv,""+endOv,""+pmStop,""+rmStop,""+spliceO,""+insFS,""+delFS,""+NumberOfLOFs,""+fullLOFs)));
		       
	    }
	   
		return new HashMap[]{SAMPLE_SUMMARY,GENE_SUMMARY};
	}
	
	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	private static void printMap(Map<String,ArrayList<String>> Map, String Outfile, String Header) {
		try {
			PrintWriter writer = new PrintWriter(Outfile, "UTF-8");
			writer.println(Header);
			
			Map<String, ArrayList<String>> map = new TreeMap<String, ArrayList<String>>(Map);
			
		    Set s = map.entrySet();
		    Iterator it = s.iterator();
		    while ( it.hasNext() ) {
		       Map.Entry entry = (Map.Entry) it.next();
		       String key = (String) entry.getKey();
		       ArrayList<String> value_list = (ArrayList<String>) entry.getValue();
		       String value = Arrays.toString(value_list.toArray());
		       value = value.replace("[", "");
		       value = value.replace("]", "");
		       value = value.replace(", ", "\t");
		       writer.println(key + "\t" + value);
		    }//while
					
		    writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}//printMap
	
	@SuppressWarnings({ "rawtypes", "unused", "unchecked" })
	private static void printMap(Map<String,ArrayList<String>> Map) {
		
		Map<String, ArrayList<String>> map = new TreeMap<String, ArrayList<String>>(Map);
		
	    Set s = map.entrySet();
	    Iterator it = s.iterator();
	    while ( it.hasNext() ) {
	       Map.Entry entry = (Map.Entry) it.next();
	       String key = (String) entry.getKey();
	       ArrayList<String> value_list = (ArrayList<String>) entry.getValue();
	       String value = Arrays.toString(value_list.toArray());
	       value = value.replace("[", "");
	       value = value.replace("]", "");
	       value = value.replace(", ", "\t");
	       System.out.println(key + "\t" + value);
	    }//while
	    System.out.println("========================");
	}//printMap
	
	
	
	private int contains(String Allele, String[] LOF_ALLELES_FIELDS){
		
//		System.out.println(Allele);
//		System.out.println(Arrays.toString(LOF_TYPE_FIELDS));
		
		for(int i = 0;i<LOF_ALLELES_FIELDS.length;i++){
			
			if(Allele.equals(LOF_ALLELES_FIELDS[i])){
				return i;
			}
		}
		return -1;
	}
	
}

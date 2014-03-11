package de.helmholtz_muenchen.ibis.ngs.denovocaller;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;

public class DeNovoCaller {

	
	private static HashMap<String, String> PED = new HashMap<String, String>();
	private static HashMap<String, String> PARENTAL_MUTATIONS = new HashMap<String, String>();
	
	
	public static BufferedDataTable[] findDeNovos(BufferedDataTable inData, String PED_FILE, String Input_Type, final ExecutionContext exec, NodeLogger LOGGER) throws Exception{
		
		/**Store the PED File**/
		readPED(PED_FILE,LOGGER);
		
		/**Check type of analysis (trio or single)**/
		if(Input_Type.equals("Single")){
			//Find the de novos
			return new BufferedDataTable[]{singleSampledeNovos(inData, exec,LOGGER)};
		}else{	//Trio
			return new BufferedDataTable[]{trioSampledeNovos(inData, exec,LOGGER)}; 
		}
	}
	
	
	
	/**
	 * Runs the deNovo Detection for single Samples
	 * @param inData that stores the input VCF files
	 * @param lOGGER 
	 */
	private static BufferedDataTable singleSampledeNovos(BufferedDataTable inData, final ExecutionContext exec, NodeLogger LOGGER){
		Iterator<DataRow> it = inData.iterator();
		LinkedList<String> ChildSamples = new LinkedList<String>();
		while(it.hasNext()){
			DataRow r = it.next();
			String VCF_FILE = r.getCell(0).toString();
			LOGGER.info("Processing File "+VCF_FILE);
			if(getSamplePhenotypeFromFile(VCF_FILE).equals("1")){	//Parental Samples
				LOGGER.info("Detected a parental sample");
				readParentalMutations(VCF_FILE);
			}else{												//Child Samples
				LOGGER.info("Detected a child sample. Storing for later processing.");
				ChildSamples.add(VCF_FILE);
			}
		}
		
		/**Processed all files, now run through the child sample files and filter for de novos**/
		
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator("DeNovoOutfiles", FileCell.TYPE).createSpec()}));
		
		Iterator<String> sampleIt = ChildSamples.iterator();
		int index = 0;	//Row Index
		while(sampleIt.hasNext()){
			String sample = sampleIt.next();
			//Filter deNovos and write to new file
			String deNovoOutfile = filterSingleDeNovos(sample);
			//Store the new files in the output table
	    	FileCell c = (FileCell) FileCellFactory.create(deNovoOutfile);
	    	cont.addRowToTable(new DefaultRow("Row"+index,c));
	    	index++;
		}
		cont.close();
		return cont.getTable();
	}
	
	/**
	 * Runs the deNovo Detection for trio Samples
	 * @param inData
	 * @param exec
	 * @param LOGGER
	 * @return
	 * @throws Exception 
	 */
	private static BufferedDataTable trioSampledeNovos(BufferedDataTable inData, final ExecutionContext exec, NodeLogger LOGGER) throws Exception{
		
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator("DeNovoOutfiles", FileCell.TYPE).createSpec()}));
		
		Iterator<DataRow> it = inData.iterator();
		int index = 0;	//Row Index
		while(it.hasNext()){
			DataRow r = it.next();
			String VCF_FILE = r.getCell(0).toString();
			LOGGER.info("Processing File "+VCF_FILE);
			String deNovoOutfile = filterTrioDeNovos(VCF_FILE);
			//Store the new files in the output table
	    	FileCell c = (FileCell) FileCellFactory.create(deNovoOutfile);
	    	cont.addRowToTable(new DefaultRow("Row"+index,c));
	    	index++;
		}
		cont.close();
		return cont.getTable();
	}
	
	
	/**
	 * Stores all mutations of a given VCF file
	 * @param VCF_FILE
	 */
	private static void readParentalMutations(String VCF_FILE){
			try(BufferedReader br = new BufferedReader(new FileReader(VCF_FILE))) {
		    	String line;
				while ((line = br.readLine()) != null) {
		            //Skips header lines
		            if(!(line.startsWith("#"))){
		            	String[] fields = line.split("\\t");
		            	//#CHROM	POS	ID	REF	ALT
		            	PARENTAL_MUTATIONS.put(fields[0]+"\t"+fields[1]+"\t"+fields[3]+"\t"+fields[4],"");		//Store each parental mutation
		            } 
		        }
		        br.close();
		    } catch (IOException e) {
		    	System.err.print("Could not find the specified VCF file !");
		        e.printStackTrace();
		    }
	}
	
	
	/**
	 * Searches for de Novo Variants and writes them to a new vcf file
	 * @param VCF_FILE
	 * @return
	 */
	private static String filterSingleDeNovos(String VCF_FILE){
		
		String deNovoVCF = replaceLast(VCF_FILE,".vcf",".deNovo.vcf");
			
		try(BufferedReader br = new BufferedReader(new FileReader(VCF_FILE))) {
			
			/**Prepare Writer and Outfile**/
			BufferedWriter writer = null;
			writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(deNovoVCF)));
			
			/**Read current VCF File**/
	    	String line;
			while ((line = br.readLine()) != null) {
	            //Writer all header lines
	            if(line.startsWith("#")){
	            	writer.write(line);
	            	writer.newLine();

	            }else{//Only write if deNovo
	            	String[] fields = line.split("\t");
	            	if(PARENTAL_MUTATIONS.containsKey(fields[0]+"\t"+fields[1]+"\t"+fields[3]+"\t"+fields[4])){
	            		// Variant is present in parental samples !!
	            	}else{
	            		// The remaining de novo variants 
	            		writer.write(line);
	            		writer.newLine();
	            	}
	            }
	        }
	        br.close();
	        writer.close();
	    } catch (IOException e) {
	    	System.err.print("Could not find the specified VCF file !");
	        e.printStackTrace();
	    }
		return deNovoVCF;
	}
	
	
	
	private static String filterTrioDeNovos(String VCF_FILE) throws Exception{
		
		int ChildIndex = -1;
		int FatherIndex = -1;
		int MotherIndex = -1;
		
		
		String deNovoVCF = VCF_FILE.replaceAll(".vcf",".deNovo.vcf");
		
		try(BufferedReader br = new BufferedReader(new FileReader(VCF_FILE))) {
			
			/**Prepare Writer and Outfile**/
			BufferedWriter writer = null;
			writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(deNovoVCF)));
			
			/**Read current VCF File**/
	    	String line;
			while ((line = br.readLine()) != null) {
	            //Writer all header lines
	            if(line.startsWith("#")){
	            	writer.write(line);
	            	writer.newLine();
	            	
	            	//Search for Header line that specifies the sample types of the multi-sample vcf file
	            	if(!(line.startsWith("##"))){	//--> Starts with # not ##
	            		String[] lineFields = line.split("\\t");
	            		for(int i = 9; i<12;i++){
	            			String tmpSampleID = lineFields[i];
	            			String sampleType = getSampleType(tmpSampleID);
	            			if(sampleType.equals("C")){
	            				ChildIndex = i;
	            			}else if(sampleType.equals("F")){
	            				FatherIndex = i;
	            			}else if(sampleType.equals("M")){
	            				MotherIndex = i;
	            			}
	            		}
	            	}

	            }else{//Only write if deNovo
	            	if(ChildIndex==-1 || MotherIndex==-1 || FatherIndex==-1){
	            	throw new Exception("Something went wrong...Probably VCF SampleIDs dont match the PED file SampleIDs!");
	            	}
	            	
	            	if(checkdeNovoTrioVariantline(line, ChildIndex, FatherIndex, MotherIndex)){
	            		// The remaining de novo variants 
	            		writer.write(line);
	            		writer.newLine();
	            	}
	            }
	        }
	        br.close();
	        writer.close();
	    } catch (IOException e) {
	    	System.err.print("Could not find the specified VCF file !");
	        e.printStackTrace();
	    }
		return deNovoVCF;
	}
	
	
	
	
	/***
	 * Reads and stores the PED file
	 * @param lOGGER 
	 * @param PED
	 */
	private static void readPED(String PED_FILE, NodeLogger LOGGER){
		
//	     Family ID
//	     Individual ID
//	     Paternal ID
//	     Maternal ID
//	     Sex (1=male; 2=female; other=unknown)
//	     Phenotype
	    
	    try(BufferedReader br = new BufferedReader(new FileReader(PED_FILE))) {
	    	LOGGER.info("Reading PED File "+PED_FILE);

	    	String line;
			while ((line = br.readLine()) != null) {
	            String[] fields = line.split("\\t");
	            String SampleID = fields[1];
	            PED.put(SampleID, line);
	        }
	        LOGGER.info("Reading PED File finished.");
	        br.close();
	    } catch (IOException e) {
	    	System.err.print("Could not find the PED file !");
	        e.printStackTrace();
	    }
	}
	
	/**
	 * Returns the stores value of the phenotype column for the provided VCF File
	 * @param VCF_FILE
	 * @return
	 */
	private static String getSamplePhenotypeFromFile(String VCF_FILE){
		String[] fields = VCF_FILE.split("/");
		String SampleID = fields[fields.length-1].split("_")[0];
		SampleID = SampleID.replaceAll(".vcf", "");			//Replace .vcf ending in case of e.g. Sample001.vcf
		String SampleType = PED.get(SampleID).split("\\t")[5];
		return SampleType;
	}
	
	
	private static String getSampleType(String SampleID) throws Exception{
		
		if(PED.get(SampleID) == null){
			try{
				throw new Exception("Something went wrong...Probably VCF SampleIDs dont match the PED file SampleIDs!");
			
			}catch (Exception e){
				System.out.println("Guessing SampleIDs...");
				SampleID = SampleID.replace("/","");
				SampleID = SampleID.split("_")[0];
				if(PED.get(SampleID) == null){
					throw new Exception("Cant guess SampleIDs...Please fix source code or PED/VCF Sample IDs !");
				}else{
					System.out.println("...worked !");
				}
			}
		}
		
		String Phenotype = PED.get(SampleID).split("\\t")[5];
		String Gender = PED.get(SampleID).split("\\t")[4];

		if(Phenotype.equals("1")){		//Parental Samples, unaffected
			if(Gender.equals("1")){
				return "F";					//Father
			}else if(Gender.equals("2")){
				return "M";					//Mother
			}else{
				return "Unknown";
			}
		}else if(Phenotype.equals("2")){
			return "C";						//Child
		}else{
			return null;
		}
	}
	
	/**
	 * Returns true if the variant is de novo
	 * @param VCF_Variantline
	 * @param ChildIndex
	 * @param FatherIndex
	 * @param MotherIndex
	 * @return
	 */
	private static boolean checkdeNovoTrioVariantline(String VCF_Variantline,int ChildIndex, int FatherIndex, int MotherIndex){
		
		String childGT = VCF_Variantline.split("\\t")[ChildIndex].split(":")[0];
		String fatherGT = VCF_Variantline.split("\\t")[FatherIndex].split(":")[0];
		String motherGT = VCF_Variantline.split("\\t")[MotherIndex].split(":")[0];
		
		if((motherGT.equals("0|0") && fatherGT.equals("0|0")) && (childGT.equals("0|1") || childGT.equals("1|0") || childGT.equals("1|1"))){
			return true;
		}else if((childGT.equals("0|1") || childGT.equals("1|1")) && fatherGT.equals("0|0")){
			return true;
		}else if((childGT.equals("1|0") || childGT.equals("1|1")) && motherGT.equals("0|0")){
			return true;
		}else if((childGT.equals("1/0") || childGT.equals("1/1") || childGT.equals("0/1") || childGT.equals("1|1") || childGT.equals("0|1") || childGT.equals("1|0")) && ( (fatherGT.equals("0|0") || fatherGT.equals("0/0")) && (motherGT.equals("0|0") || motherGT.equals("0/0")))){
			return true;
		}else{
			return false;
		}	
	}
	
    public static String replaceLast(String original, String target, String replacement) {
        int index = original.lastIndexOf(target);
       
        if(index == -1) {
                return original;
        }
       
        return original.substring(0, index) + replacement + original.substring(index + target.length());
}
	
}

package de.helmholtz_muenchen.ibis.ngs.vcffilter;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

import org.knime.core.node.NodeLogger;

public class VATFilter {
	
	//SO_term -> VAT term
	HashMap<String,String> SOtoVAT;
	
	HashSet<String> lookUpTerms;
	
	StringBuilder header;
	StringBuilder filteredContent;
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(VCFFilterNodeModel.class);
	
	public VATFilter() {
		SOtoVAT = new HashMap<>();
		SOtoVAT.put("splice_donor_variant", "spliceOverlap");
		SOtoVAT.put("splice_acceptor_variant", "spliceOverlap");
		SOtoVAT.put("stop_gained", "prematureStop");
		SOtoVAT.put("frameshift_variant", "insertionFS");
		SOtoVAT.put("frameshift_variant", "deletionFS");
		SOtoVAT.put("stop_lost","removedStop");
		SOtoVAT.put("initiator_codon_variant","startOverlap");
		SOtoVAT.put("inframe_insertion", "insertionNFS");
		SOtoVAT.put("inframe_deletion", "deletionNFS");
		SOtoVAT.put("missense_variant", "nonsynonymous");
	}

	public void filter(String infile, String outfile, HashSet<String> SO_terms) {
		
		//fill lookUpTerms
		lookUpTerms = new HashSet<String>();
		for(String t: SO_terms) {
			if(SOtoVAT.containsKey(t)) {
				lookUpTerms.add(SOtoVAT.get(t));
			} else {
				LOGGER.warn("No VAT term found for SO term:"+t);
			}
		}
		
		readFile(infile);
		writeFile(outfile);
	}
	
	private void readFile(String infile) {
		
		//read file
		header = new StringBuilder();
		filteredContent = new StringBuilder();
		int infoPos = -1;
		String [] fields;
		try {
			FileInputStream inputStream = null;
			Scanner sc = null;
			try {
			    inputStream = new FileInputStream(infile);
			    sc = new Scanner(inputStream, "UTF-8");
			    while (sc.hasNextLine()) {
			        String line = sc.nextLine();
			        if(line.startsWith("#"))  {
			        	header.append(line + System.getProperty("line.separator"));
			        	if(line.startsWith("#CHROM")) {
			        		fields = line.split("\t");
			        		for(int i = 0; i < fields.length; i++) {
			        			if(fields[i].equals("INFO")) {
			        				infoPos = i;
			        				break;
			        			}
			        		}
			        	}
			        } else {
			        	//do filtering
			        	String filteredLine = "";
			        	fields = line.split("\t");
			        	
			        	String [] info_fields = fields[infoPos].split(";");
			        	String [] allele_annotations = null;
			        	int VApos = -1;
			        	
			        	for(int i = 0; i < info_fields.length; i++) {
			        		String info = info_fields[i];
			    			if(info.startsWith("VA")) {
			    				info = info.replaceFirst("VA=","");
			    				allele_annotations = info.split(",");
			    				VApos = i;
			    				break;
			    			}
			    		}
			    		
			    		if(allele_annotations == null) continue; //nothing was annotated
			    		
			    		ArrayList<String> passed_annotations = new ArrayList<>();
			    		for(String a: allele_annotations) {
			    			String consequence = a.split(":")[4];
			    			if(lookUpTerms.contains(consequence)) {
			    				passed_annotations.add(a);
			    			}
			    		}
			    		
			    		//fill filteredLine
			    		if(passed_annotations.size()>0) {
			    			//append fields before INFO field
			    			for(int i = 0; i < infoPos; i++) {
			    				filteredLine += fields[i] + "\t";
			    			}
			    			
			    			//append info fields before VA
			    			for(int i = 0; i < VApos; i++) {
			    				filteredLine += info_fields[i] + ";";
			    			}
			    			
			    			//append passed annotations
			    			filteredLine += "VA="+passed_annotations.get(0);
			    			for(int j = 1; j<passed_annotations.size(); j++) {
			    				filteredLine += "," + passed_annotations.get(j);
			    			}
			    			
			    			//append info fields after VA
			    			for(int i = VApos+1; i < info_fields.length; i++) {
			    				filteredLine += ";" + info_fields[i];
			    			}
			    			
			    			//append fields after INFO field
			    			for(int k = infoPos+1; k< fields.length; k++) {
			    				filteredLine += "\t" + fields[k];
			    			}
			    			filteredLine += System.getProperty("line.separator");
			    			filteredContent.append(filteredLine);
			    		}
			    		
			        }
			    }
			    // note that Scanner suppresses exceptions
			    if (sc.ioException() != null) {
			        throw sc.ioException();
			    }
			} finally {
			    if (inputStream != null) {
			        inputStream.close();
			    }
			    if (sc != null) {
			        sc.close();
			    }
			}
			
		} catch (IOException e) {
			LOGGER.error("VCF file couldn't be read!");
			e.printStackTrace();
		}
	}
	
	private void writeFile(String outfile) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
			bw.write(header.toString());
			bw.write(filteredContent.toString());
			bw.close();
		} catch (IOException e) {
			LOGGER.error("Filtered VCF file couldn't be written!");
			e.printStackTrace();
		}
	}
	
//	public static void main (String [] args) {
//		String infile = "/home/ibis/tim.jeske/KORAdata/analysis_ready/20150319_vat_k/analysis_ready.diabetes.filtered.haplotypecaller.vat.vcf";
//		String outfile = infile.replace("vcf", "filtered.vcf");
//		HashSet<String> myTerms = new HashSet<>();
//		myTerms.add("stop_gained");
//		new VATFilter().filter(infile,outfile, myTerms);
//	}
}

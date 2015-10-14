package de.helmholtz_muenchen.ibis.utils.ngs;

import java.util.HashMap;
import java.util.HashSet;

public class VEPAnnotationParser {
	
	public static final String INFO_ID = "CSQ";
	
	int allele_num_index = -1;
	int gene_id_index = -1;
	
	/**
	 * 
	 * @param vep_header description of the INFO field for VEP annotations, can be retrieved by getInfoHeader(VEPAnnotationParser.INFO_ID) for a VCFFile 
	 * @throws Exception is thrown if vep_header does not start with the ID of VEP annotations
	 */
	public VEPAnnotationParser(String vep_header) throws Exception {
		
		if(!vep_header.startsWith("ID="+INFO_ID)) {
			throw new IllegalArgumentException("Invalid VEP header! It does not start with ID="+INFO_ID+"!");
		}
		
		String [] vep_fields = vep_header.split("\\|");
    	for(int i = 0; i < vep_fields.length; i++) {
    		if(vep_fields[i].equals("ALLELE_NUM")) {
    			this.allele_num_index = i;
    		} else if(vep_fields[i].equals("Gene")) {
    			this.gene_id_index = i;
    		}
    	}	
	}
	
	/**
	 * 
	 * @param anno VEP annotation string in the INFO field of a variant, can be retrieved by getInfoField(VEPAnnotationParser.INFO_ID) for a VCFVariant
	 * @return annotated alleles (ALLELE_NUM) for each gene affected by a variant
	 */
	public HashMap<String, HashSet<Integer>> getGene2AlleleNum (String anno) {
		if(allele_num_index == -1 || gene_id_index == -1) {
			throw new IllegalArgumentException("ALLELE_NUM and/or gene have not been annotated by VEP!");
		}
		
		HashMap<String, HashSet<Integer>> gene2allele_num = new HashMap<>();
		int allele_num;
		String gene = "";
		
		//iterate over all transcript annotations
		for(String a : anno.split(",")) {
			gene = a.split("\\|")[gene_id_index];
			allele_num = Integer.parseInt(a.split("\\|")[allele_num_index]);
			if(gene2allele_num.containsKey(gene)) {
				gene2allele_num.get(gene).add(allele_num);
			} else {
				HashSet<Integer> tmp = new HashSet<>();
				tmp.add(allele_num);
				gene2allele_num.put(gene,tmp);
			}
		}
		return gene2allele_num;
	}
	

}

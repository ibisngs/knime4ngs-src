package de.helmholtz_muenchen.ibis.utils.ngs;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import de.helmholtz_muenchen.ibis.ngs.lofsummary.Annotation;

public class VEPAnnotationParser implements AnnotationParser {
	
	public final static String ANN_ID = "CSQ";
	
	int allele_id_index = -1;
	int gene_id_index = -1;
	int transcript_id_index = -1;
	int consequence_index = -1;
	int symbol_index = -1;
	
	/**
	 * 
	 * @param vep_header description of the INFO field for VEP annotations, can be retrieved by getInfoHeader(VEPAnnotationParser.INFO_ID) for a VCFFile 
	 * @throws Exception is thrown if vep_header does not start with the ID of VEP annotations
	 */
	public VEPAnnotationParser(String vep_header) throws Exception {
		
		if(!vep_header.startsWith("ID="+getAnnId())) {
			throw new IllegalArgumentException("Invalid VEP header! It does not start with ID="+getAnnId()+"!");
		}
		
		String [] vep_fields = vep_header.split("\\|");
		String part;
    	for(int i = 0; i < vep_fields.length; i++) {
    		part = vep_fields[i];
    		if(part.equals("ALLELE_NUM")) {
    			this.allele_id_index = i;
    		} else if(part.equals("Gene")) {
    			this.gene_id_index = i;
    		} else if(part.equals("Feature")) {
    			this.transcript_id_index = i;
    		} else if(part.contains("Consequence")) {
    			this.consequence_index = i;
    		} else if(part.equals("SYMBOL")) {
    			this.symbol_index = i;
    		}
    	}	
	}
	
	/**
	 * 
	 * @param anno VEP annotation string in the INFO field of a variant, can be retrieved by getInfoField(VEPAnnotationParser.INFO_ID) for a VCFVariant
	 * @return annotated alleles (ALLELE_NUM) for each gene affected by a variant
	 */
	public HashMap<String, HashSet<Integer>> getGene2AlleleIds (String anno, boolean use_id) {
		if(allele_id_index == -1 || gene_id_index == -1) {
			throw new IllegalArgumentException("ALLELE_NUM and/or gene have not been annotated by VEP!");
		}
		
		HashMap<String, HashSet<Integer>> gene2allele_num = new HashMap<>();
		int allele_num;
		String gene = "";
		
		//iterate over all transcript annotations
		for(String a : anno.split(",")) {
			if(use_id) {
				gene = a.split("\\|")[gene_id_index];
			} else {
				gene = a.split("\\|")[symbol_index];
			}
			
			allele_num = Integer.parseInt(a.split("\\|")[allele_id_index]);
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

	@Override
	public String getAnnId() {
		return ANN_ID;
	}

	@Override
	public HashMap<String, HashSet<String>> getGene2TranscriptIds(int allele_id, String anno) {
		if(allele_id_index == -1 || gene_id_index == -1 || transcript_id_index == -1) {
			throw new IllegalArgumentException("ALLELE_NUM and/or gene and/or feature have not been annotated by VEP!");
		}
		HashMap<String, HashSet<String>> res = new HashMap<>();
		int allele_num;
		String gene_id, transcript_id;
		String [] fields;
		for(String a: anno.split(",")) {
			fields = a.split("\\|");
			allele_num = Integer.parseInt(fields[allele_id_index]);
			if(allele_num == allele_id) {
				gene_id = fields[gene_id_index];
				transcript_id = fields[transcript_id_index];
				if(res.containsKey(gene_id)) {
					res.get(gene_id).add(transcript_id);
				} else {
					HashSet<String> tmp = new HashSet<>();
					tmp.add(transcript_id);
					res.put(gene_id, tmp);
				}
			}
		}
		return res;
	}
	
	@Override
	public HashSet<String> getAnnotatedAlleleIds(String anno) {
		if(allele_id_index == -1) {
			throw new IllegalArgumentException("ALLELE_NUM has not been annotated!");
		}
		HashSet<String> result = new HashSet<>();
		String allele_id;
		for(String a : anno.split(",")) {
			allele_id = a.split("\\|")[allele_id_index];
			result.add(allele_id);
		}
		return result;
	}

	@Override
	public LinkedList<Annotation> getAnnotations(int allele_id, String anno) {
		if(allele_id_index == -1 || gene_id_index == -1 || transcript_id_index == -1 || consequence_index == -1) {
			throw new IllegalArgumentException("ALLELE_NUM and/or gene and/or feature and/or consequence have not been annotated by VEP!");
		}
		LinkedList<Annotation> res = new LinkedList<>();
		int allele_num;
		String gene_id, transcript_id, consequence;
		String [] fields;
		for(String a: anno.split(",")) {
			fields = a.split("\\|");
			allele_num = Integer.parseInt(fields[allele_id_index]);
			if(allele_num == allele_id) {
				gene_id = fields[gene_id_index];
				transcript_id = fields[transcript_id_index];
				consequence = fields[consequence_index];
				res.add(new Annotation(gene_id, transcript_id, consequence));
			}
		}
		return res;
	}
}

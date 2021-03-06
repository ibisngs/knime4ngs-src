/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.utils.ngs;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import de.helmholtz_muenchen.ibis.ngs.vepsummary.Annotation;

public class VEPAnnotationParser implements AnnotationParser {
	
	public final static String ANN_ID = "CSQ";
	
	int allele_id_index = -1;
	int gene_id_index = -1;
	int transcript_id_index = -1;
	int consequence_index = -1;
	int symbol_index = -1;
	int exac_af_NFE_index = -1;
	
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
    		} else if(part.equals("ExAC_AF_NFE")) {
    			this.exac_af_NFE_index = i;
    		}
    	}	
	}
	
	/**
	 * 
	 * @param anno VEP annotation string in the INFO field of a variant, can be retrieved by getInfoField(VEPAnnotationParser.INFO_ID) for a VCFVariant
	 * @return annotated alleles (ALLELE_NUM) for each entity affected by a variant
	 */
	public HashMap<String, HashSet<Integer>> getEntity2AlleleIds (String anno, BioEntity e) {
		
		HashMap<String, HashSet<Integer>> entity2allele_num = new HashMap<>();
		int allele_num;
		String ent = "";
		
		//iterate over all transcript annotations
		for(String a : anno.split(",")) {
			switch(e) {
				case GENE_ID:
					ent = a.split("\\|")[gene_id_index];
					break;
				case GENE_SYMBOL:
					ent = a.split("\\|")[symbol_index];
					break;
				case TRANSCRIPT_ID:
					ent = a.split("\\|")[transcript_id_index];
					break;
			}
			
			if(ent.equals("") || ent == null) {
				continue;
			}
			
			allele_num = Integer.parseInt(a.split("\\|")[allele_id_index]);
			if(entity2allele_num.containsKey(ent)) {
				entity2allele_num.get(ent).add(allele_num);
			} else {
				HashSet<Integer> tmp = new HashSet<>();
				tmp.add(allele_num);
				entity2allele_num.put(ent,tmp);
			}
		}
		return entity2allele_num;
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
				if(gene_id.equals("") || gene_id == null) {
					continue;
				}
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
	
	public double getExACAF(int allele_id, String anno) {
		if(allele_id_index == -1 || exac_af_NFE_index == -1) {
			return 0.0;
		}
		String [] fields;
		int allele_num;
		for(String a: anno.split(",")) {
			fields = a.split("\\|");
			allele_num = Integer.parseInt(fields[allele_id_index]);
			if(allele_num == allele_id) {
				try {
					return Double.parseDouble(fields[exac_af_NFE_index]);
				} catch (Exception e) {
					return 0.0;
				}
			}
		}
		return 0.0;
	}
}

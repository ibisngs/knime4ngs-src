package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import de.helmholtz_muenchen.ibis.ngs.caseControlAnalyzer.ContingencyTable;
import de.helmholtz_muenchen.ibis.utils.ngs.AnnotationParser;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFFile;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFVariant;

public class RegionSummary {
	
	ArrayList<String> sample_ids;
	HashMap<String, Set<String>> unaffected_samples = new HashMap<>();
	HashMap<String, Set<String>> affected_samples = new HashMap<>();
	HashMap<String, Set<String>> enrich_unaffected_samples = new HashMap<>();
	HashMap<String, Set<String>> enrich_affected_samples = new HashMap<>();
	
	/**
	 * Computes a set of unaffected and affected samples for each annotated gene in the given vcf file.
	 * @param vcf
	 * @param ap
	 */
	public RegionSummary(VCFFile vcf, AnnotationParser ap, boolean use_id) {
		unaffected_samples = new HashMap<>();
		affected_samples = new HashMap<>();
		sample_ids = vcf.getSampleIds();
		
		HashMap<String, HashSet<Integer>> gene2allele_ids;
		HashSet<Integer> alt_alleles;
		HashSet<String> aff_samp, unaff_samp;
		VCFVariant var;
		String anno;
		
		while(vcf.hasNext()) {
			var = vcf.next();
			anno = var.getInfoField(ap.getAnnId());
			gene2allele_ids = ap.getGene2AlleleIds(anno, use_id); //use gene symbol
			
			for(String gene: gene2allele_ids.keySet()) {
				alt_alleles = gene2allele_ids.get(gene);
				aff_samp = var.getAffectedSamples(alt_alleles);
				unaff_samp = var.getUnaffectedSamples(alt_alleles);
				if(unaffected_samples.containsKey(gene)) {
					//unaffected samples contains all samples that carry only non-LOF alleles in a gene
					unaffected_samples.get(gene).retainAll(unaff_samp);
				} else {
					unaffected_samples.put(gene,unaff_samp);
				}
				
				if(affected_samples.containsKey(gene)) {
					affected_samples.get(gene).addAll(aff_samp);
				} else {
					affected_samples.put(gene, aff_samp);
				}
			}
		}
	}
	
	private HashMap<String, Double> getFrequencies(HashMap<String, Set<String>> aff, HashMap<String, Set<String>> un) {
		HashMap<String, Double> result = new HashMap<>();
		
		HashSet<String> called_samples;
		
		for(String gene: aff.keySet()) {
			called_samples = new HashSet<>();
			called_samples.addAll(un.get(gene));
			called_samples.addAll(aff.get(gene));
			result.put(gene, (double)aff.get(gene).size()/(double)called_samples.size());
		}
	
		return result;
	}
	
	public HashMap<String, ContingencyTable> getTables(HashMap<String, Boolean> sample_id2affected, HashMap<String, Set<String>> aff, HashMap<String, Set<String>> un) {
		HashMap<String, ContingencyTable> result = new HashMap<>();
		
		int aff_case, aff_ctrl, un_case, un_ctrl;
		for(String gene: aff.keySet()) {
			
			aff_case = 0;
			aff_ctrl = 0;
			un_case = 0;
			un_ctrl = 0;
			
			for(String s: aff.get(gene)) {
				if(sample_id2affected.containsKey(s) && sample_id2affected.get(s)) {
					aff_case++;
				} else {
					aff_ctrl++;
				}
			}
			
			for(String s: un.get(gene)) {
				if(sample_id2affected.containsKey(s) && sample_id2affected.get(s)) {
					un_case++;
				} else {
					un_ctrl++;
				}
			}
			
			if(aff_case>0 || aff_ctrl>0) {
				result.put(gene, new ContingencyTable(aff_case, un_case, aff_ctrl, un_ctrl));
			}
		}
		return result;
	}
	
	
	/**
	 * Computes #affected_samples/#called_samples for each gene.
	 * @return a map gene -> frequency
	 */
	public HashMap<String, Double> getFrequencies() {
		return getFrequencies(affected_samples, unaffected_samples);
	}
	
	/**
	 * Computes #affected_samples/#called_samples for each set definition.
	 * @return a map gene -> frequency
	 */
	public HashMap<String, Double> getSetFrequencies() {
		return getFrequencies(enrich_affected_samples, enrich_unaffected_samples);
	}
	
	public HashMap<String, ContingencyTable> getTables(HashMap<String, Boolean> sample_id2affected) {
		return getTables(sample_id2affected, affected_samples, unaffected_samples);
	}
	
	public HashMap<String, ContingencyTable> getSetTables(HashMap<String, Boolean> sample_id2affected) {
		return getTables(sample_id2affected, enrich_affected_samples, enrich_unaffected_samples);
	}
	
	public void groupBy(HashMap<String, HashSet<String>> id2set) {
		HashSet<String> genes;
		
		for(String set:id2set.keySet()) {
			genes = id2set.get(set);
			
			Set<String> tmp_enrich_un = new HashSet<>();
			for(String s: sample_ids) {
				tmp_enrich_un.add(s);
			}
			Set<String> tmp_enrich_aff = new HashSet<>();
			for(String g: genes) {
				if(unaffected_samples.containsKey(g)) { 
					tmp_enrich_un.retainAll(unaffected_samples.get(g));
					tmp_enrich_aff.addAll(affected_samples.get(g));
				}
			}
			enrich_unaffected_samples.put(set, tmp_enrich_un);
			enrich_affected_samples.put(set, tmp_enrich_aff);
		}
	}
}

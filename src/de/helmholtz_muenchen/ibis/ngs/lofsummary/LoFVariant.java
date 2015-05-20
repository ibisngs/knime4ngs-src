package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.ArrayList;
import java.util.HashMap;

public class LoFVariant {
	
	String chr, pos, ref_allele, alt_allele, rsId, af;
	int observed_homo, observed_hetero;
	HashMap<String, LoFGene> gene_id2lofgene_fields;
	ArrayList<String> genotypes;
	
	public LoFVariant (String chr, String pos, String ref_allele, String alt_allele, String rsId, String af, int observed_homo, int observed_hetero, HashMap<String, LoFGene> gene_id2lofgene_fields, ArrayList<String> genotypes) {
		this.chr = chr.replaceAll("[chroCHRO]","");
		this.pos = pos;
		this.ref_allele = ref_allele;
		this.alt_allele = alt_allele;
		this.af = af;
		this.rsId = rsId;
		this.observed_homo = observed_homo;
		this.observed_hetero = observed_hetero;
		this.gene_id2lofgene_fields = gene_id2lofgene_fields;
		this.genotypes = genotypes;
	}

	public String getChr() {
		return chr;
	}

	public String getPos() {
		return pos;
	}

	public String getRef_allele() {
		return ref_allele;
	}

	public String getAlt_allele() {
		return alt_allele;
	}

	public String getRsId() {
		return rsId;
	}
	
	public String getAF() {
		return af;
	}

	public int getObserved_homo() {
		return observed_homo;
	}

	public int getObserved_hetero() {
		return observed_hetero;
	}

	public HashMap<String, LoFGene> getGene_id2lofgene_fields() {
		return gene_id2lofgene_fields;
	}
	
	public ArrayList<String> getGenotypes() {
		return genotypes;
	}

}

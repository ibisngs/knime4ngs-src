package de.helmholtz_muenchen.ibis.ngs.lofstatistics;

import java.util.ArrayList;
import java.util.HashMap;

public class LoFVariant {
	
	String chr, pos, ref_allele, alt_allele, rsId;
	int observed_homo, observed_hetero;
	HashMap<String, LoFGene> gene_id2lofgene_fields;
	ArrayList<String> genotypes;
	
	public LoFVariant (String chr, String pos, String ref_allele, String alt_allele, String rsId, int observed_homo, int observed_hetero, HashMap<String, LoFGene> gene_id2lofgene_fields, ArrayList<String> genotypes) {
		this.chr = chr.replaceAll("[chroCHRO]","");
		this.pos = pos;
		this.ref_allele = ref_allele;
		this.alt_allele = alt_allele;
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

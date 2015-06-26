package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.ArrayList;

public class SampleInfo {
	
	String fam_id, pat_id, mat_id, sex, phenotype;
	int fullLOFs, partLOFs;
	ArrayList<String> complete_LOF_genes;
	ArrayList<String> part_LOF_genes;
	ArrayList<String> hom_LOF_genes;
	
	public SampleInfo() {
		this.fam_id = "unknown";
		this.pat_id = "unknown";
		this.mat_id = "unknown";
		this.sex = "unknown";
		this.phenotype = "unknown";
		fullLOFs = 0;
		partLOFs = 0;
		complete_LOF_genes = new ArrayList<>();
		part_LOF_genes = new ArrayList<>();
		hom_LOF_genes = new ArrayList<>();
	}
	
	public SampleInfo(String fam_id, String pat_id, String mat_id, String sex, String phenotype) {
		this.fam_id = fam_id;
		this.pat_id = pat_id;
		this.mat_id = mat_id;
		this.sex = sex;
		this.phenotype = phenotype;
		fullLOFs = 0;
		partLOFs = 0;
		complete_LOF_genes = new ArrayList<>();
		part_LOF_genes = new ArrayList<>();
		hom_LOF_genes = new ArrayList<>();
	}
	
	public String getPatId() {
		return pat_id;
	}
	
	public String getMatId() {
		return mat_id;
	}
	
	public int getFullLOFs() {
		return fullLOFs;
	}

	public int getPartLOFs() {
		return partLOFs;
	}

	public void incrementFullLOFs () {
		fullLOFs++;
	}
	
	public void incrementPartLOFs () {
		partLOFs++;
	}
	
	public void addCompleteLoFGene(String gene_id) {
		complete_LOF_genes.add(gene_id);
	}
	
	public void addPartLoFGene(String gene_id) {
		part_LOF_genes.add(gene_id);
	}
	
	public void addHomLoFGene(String gene_id) {
		hom_LOF_genes.add(gene_id);
	}
	
	public ArrayList<String> getComplete_LOF_genes() {
		return complete_LOF_genes;
	}

	public ArrayList<String> getPart_LOF_genes() {
		return part_LOF_genes;
	}
	
	public ArrayList<String> getHom_LOF_genes() {
		return this.hom_LOF_genes;
	}
}

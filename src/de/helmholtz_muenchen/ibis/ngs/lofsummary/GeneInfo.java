package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.HashSet;


public class GeneInfo {
	
	String contig;
	String symbol;
	
	HashSet<String> unaffected_samples;
	HashSet<String> affected_samples;
	HashSet<String> hom_samples;
	HashSet<String> ko_samples;
	
	int fullLoFs, partLoFs, aff_case, aff_ctrl, un_case, un_ctrl;
	
	public GeneInfo() {
		fullLoFs = 0;
		partLoFs = 0;
		aff_case = 0;
		aff_ctrl = 0; 
		un_case = 0; 
		un_ctrl = 0;
		contig = "";
		affected_samples = new HashSet<>();
		ko_samples = new HashSet<>();
		hom_samples = new HashSet<>();
		unaffected_samples = new HashSet<>();
	}
	
	public int getAff_case() {
		return aff_case;
	}

	public int getAff_ctrl() {
		return aff_ctrl;
	}

	public int getUn_case() {
		return un_case;
	}

	public int getUn_ctrl() {
		return un_ctrl;
	}
	
	public void setTable(int aff_case, int aff_ctrl, int un_case, int un_ctrl) {
		this.aff_case = aff_case;
		this.aff_ctrl = aff_ctrl;
		this.un_case = un_case;
		this.un_ctrl = un_ctrl;
	}

	public void setContig(String contig) {
		this.contig = contig;
	}
	
	public String getContig() {
		return this.contig;
	}
	
	public void setSymbol(String name) {
		this.symbol = name;
	}
	
	public String getSymbol() {
		return this.symbol;
	}
	
	public void incrementFullLoFs () {
		fullLoFs++;
	}
	
	public int getFullLoFs () {
		return this.fullLoFs;
	}
	
	public void incrementPartLoFs () {
		partLoFs++;
	}
	
	public int getPartLoFs() {
		return this.partLoFs;
	}
	

	public void addUnaffectedSample(String sample_id) {
		this.unaffected_samples.add(sample_id);
	}
	
	public HashSet<String> getUnaffectedSamples() {
		return this.unaffected_samples;
	}
	
	public void addAffectedSample(String sample_id) {
		this.affected_samples.add(sample_id);
	}
	
	public HashSet<String> getAffectedSamples () {
		return this.affected_samples;
	}
	
	public void addKOSample(String sample_id) {
		this.ko_samples.add(sample_id);
	}
	
	public HashSet<String> getKOSamples() {
		return this.ko_samples;
	}
	
	public HashSet<String> getHom_samples() {
		return hom_samples;
	}

	public void addHomSample(String sample_id) {
		this.hom_samples.add(sample_id);
	}
}

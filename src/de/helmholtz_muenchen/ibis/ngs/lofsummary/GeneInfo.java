package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.HashMap;
import java.util.HashSet;


public class GeneInfo implements Comparable<GeneInfo> {
	
	String contig;
	String symbol;

	HashMap<String, Double> pos2af_prob;
	
	HashSet<String> unaffected_samples;
	HashSet<String> affected_samples;
	HashSet<String> hom_samples;
	HashSet<String> ko_samples;
	double p_lof_aff, p_val_case_vs_control, p_val_vs_bg, p_val_fisher;
	
	int fullLoFs, partLoFs, aff_case, aff_ctrl, un_case, un_ctrl;
	
	public GeneInfo() {
		fullLoFs = 0;
		partLoFs = 0;
		aff_case = 0;
		aff_ctrl = 0; 
		un_case = 0; 
		un_ctrl = 0;
		p_val_fisher = -1.0;
		p_val_case_vs_control = -1.0;
		p_val_vs_bg = -1.0;
		p_lof_aff = -1.0;
		contig = "";
		pos2af_prob = new HashMap<>();
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
	
	public void calcFisher(int aff_case, int aff_ctrl, int un_case, int un_ctrl) {
		this.aff_case = aff_case;
		this.aff_ctrl = aff_ctrl;
		this.un_case = un_case;
		this.un_ctrl = un_ctrl;
		p_val_fisher = FisherExact.oneTailedGreater(aff_case, aff_ctrl, un_case, un_ctrl);
	}
	
	public double getP_val_Fisher() {
		return this.p_val_fisher;
	}
	
	public double getP_lof_aff() {
		return p_lof_aff;
	}

	public void setP_lof_aff(double p_lof_aff) {
		this.p_lof_aff = p_lof_aff;
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
	
	public double getProbUnaffected() {
		double result = 1.0;
		for(String pos: pos2af_prob.keySet()) {
			result = result * Math.pow(1.0 - pos2af_prob.get(pos),2);
		}
		return result;
	}
	
	public void addProb(String pos, double prob) {
		if(pos2af_prob.containsKey(pos)) {
			pos2af_prob.put(pos, prob+pos2af_prob.get(pos));
		} else {
			pos2af_prob.put(pos, prob);
		}
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

	@Override
	public int compareTo(GeneInfo gi) {
		return Double.compare(p_val_fisher, gi.getP_val_Fisher());
	}

	public double getP_val_case_vs_control() {
		return p_val_case_vs_control;
	}

	public void setP_val_case_vs_control(double p_val_case_vs_control) {
		this.p_val_case_vs_control = p_val_case_vs_control;
	}

	public double getP_val_vs_exac() {
		return p_val_vs_bg;
	}

	public void setP_val_vs_bg(double p_val_vs_exac) {
		this.p_val_vs_bg = p_val_vs_exac;
	}
}

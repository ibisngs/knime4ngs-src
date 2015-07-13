package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

public class GeneInfo {
	
	String contig;
	String symbol;
	int fullLoFs;
	int partLoFs;
	HashMap<String, Double> pos2af_prob;
	HashSet<String> affected_samples;
	HashSet<String> ko_samples;
	HashSet<String> unaffected_samples;
	
	public GeneInfo() {
		fullLoFs = 0;
		partLoFs = 0;
		contig = "";
		pos2af_prob = new HashMap<>();
		affected_samples = new HashSet<>();
		ko_samples = new HashSet<>();
		unaffected_samples = new HashSet<>();
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
	
	public double getStdDev() {
		SummaryStatistics sum = new SummaryStatistics();
		for(double value: pos2af_prob.values()) {
			sum.addValue(value);
		}
		return sum.getStandardDeviation();
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
}

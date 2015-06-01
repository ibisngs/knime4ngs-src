package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.HashMap;

public class GeneInfo {
	
	String symbol;
	int fullLoFs;
	int partLoFs;
	HashMap<String, Double> pos2af_prob;
	
	public GeneInfo() {
		fullLoFs = 0;
		partLoFs = 0;
		pos2af_prob = new HashMap<>();
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
		pos2af_prob.put(pos, prob);
	}

}

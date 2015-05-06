package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.ArrayList;

public class LoFGene {

	protected ArrayList<String> consequences;
	protected ArrayList<String> transcripts;
	protected ArrayList<ArrayList<String>> additional_infos;

	public LoFGene(int n) {
		consequences = new ArrayList<>();
		transcripts = new ArrayList<>();
		additional_infos = new ArrayList<ArrayList<String>>();
		for(int i =0; i< n; i++) {
			additional_infos.add(new ArrayList<String>());
		}
	}

	public void addConsequence(String c) {
		this.consequences.add(c);
	}

	public void addTranscript(String t) {
		this.transcripts.add(t);
	}
	
	public void addInfo(String info, int index) {
		this.additional_infos.get(index).add(info);
	}
	
	public ArrayList<String> getInfo(int index) {
		return this.additional_infos.get(index);
	}

	public int getNrOfTranscripts() {
		return transcripts.size();
	}

	public ArrayList<String> getConsequences() {
		return consequences;
	}

	public ArrayList<String> getTranscripts() {
		return transcripts;
	}
}
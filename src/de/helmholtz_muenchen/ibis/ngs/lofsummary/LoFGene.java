package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.ArrayList;
import java.util.HashSet;

public class LoFGene {

	protected HashSet<String> consequences;
	protected HashSet<String> transcripts;
	protected ArrayList<HashSet<String>> additional_infos;

	public LoFGene(int n) {
		consequences = new HashSet<String>();
		transcripts = new HashSet<String>();
		additional_infos = new ArrayList<HashSet<String>>();
		for(int i =0; i< n; i++) {
			additional_infos.add(new HashSet<String>());
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
	
	public HashSet<String> getInfo(int index) {
		return this.additional_infos.get(index);
	}

	public int getNrOfTranscripts() {
		return transcripts.size();
	}

	public HashSet<String> getConsequences() {
		return consequences;
	}

	public HashSet<String> getTranscripts() {
		return transcripts;
	}
}
package de.helmholtz_muenchen.ibis.ngs.lofstatistics;

import java.util.ArrayList;

public class LoFGene {
	
	ArrayList<String> consequences, transcripts, confidences, filter, infos, flags;
	
	public LoFGene() {
		consequences = new ArrayList<String>();
		transcripts = new ArrayList<String>();
		confidences = new ArrayList<String>();
		filter = new ArrayList<String>();
		infos = new ArrayList<String>();
		flags = new ArrayList<String>();
	}
	
	public void addConsequence(String c) {
		this.consequences.add(c);
	}
	
	public void addTranscript(String t) {
		this.transcripts.add(t);
	}
	
	public void addConfidences(String c) {
		this.confidences.add(c);
	}
	
	public void addFilter(String f) {
		this.filter.add(f);
	}
	
	public void addInfo(String i) {
		this.infos.add(i);
	}
	
	public void addFlag(String f) {
		this.flags.add(f);
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

	public ArrayList<String> getConfidences() {
		return confidences;
	}

	public ArrayList<String> getFilter() {
		return filter;
	}

	public ArrayList<String> getInfos() {
		return infos;
	}

	public ArrayList<String> getFlags() {
		return flags;
	}
	
}

package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.HashSet;

public class GeneSampleInfo {
	
	HashSet<String> transcripts_hom;
	HashSet<String> transcripts_het;
	HashSet<String> transcripts_pat;
	HashSet<String> transcripts_mat;
	boolean unaffected;
	int het_trans, hom_trans;
	
	public GeneSampleInfo() {
		unaffected = true;
		transcripts_hom = new HashSet<>();
		transcripts_het = new HashSet<>();
		transcripts_pat = new HashSet<>();
		transcripts_mat = new HashSet<>();
	}
	
	public void setAffected() {
		this.unaffected = false;
	}
	
	public boolean isUnaffected() {
		return unaffected;
	}
	
	public void addHomTranscripts(HashSet<String> transcript_ids) {
		this.transcripts_hom.addAll(transcript_ids);
	}

	public void addHetTranscripts(HashSet<String> transcript_ids) {
		this.transcripts_het.addAll(transcript_ids);
	}
	
	public void addPatTranscripts(HashSet<String> transcript_ids) {
		this.transcripts_pat.addAll(transcript_ids);
	}
	
	public void addMatTranscripts(HashSet<String> transcript_ids) {
		this.transcripts_mat.addAll(transcript_ids);
	}
	
	public void calculate() {
//		HashSet<String> hom = new HashSet<>();
//		hom.addAll(transcripts_hom);
		HashSet<String> intersection = new HashSet<String>(transcripts_mat);
		intersection.retainAll(transcripts_pat);
//		hom.addAll(intersection);
		transcripts_hom.addAll(intersection);
//		this.hom_trans = hom.size();
		this.hom_trans = transcripts_hom.size();
		
//		HashSet<String> het = new HashSet<>();
//		het.addAll(transcripts_het);
//		het.addAll(transcripts_mat);
		transcripts_het.addAll(transcripts_mat);
//		het.addAll(transcripts_pat);
		transcripts_het.addAll(transcripts_pat);
//		het.removeAll(hom);
		transcripts_het.removeAll(transcripts_hom);
		this.het_trans = transcripts_het.size();
		
//		transcripts_hom.clear();
//		transcripts_het.clear();
//		transcripts_pat.clear();
//		transcripts_mat.clear();
	}
	
	public int getHomTrans() {
		return hom_trans;
	}
	
	public int getHetTrans() {
		return het_trans;
	}

	public HashSet<String> getTranscripts_hom() {
		return transcripts_hom;
	}

	public HashSet<String> getTranscripts_het() {
		return transcripts_het;
	}

}

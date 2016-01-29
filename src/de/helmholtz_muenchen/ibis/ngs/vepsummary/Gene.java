package de.helmholtz_muenchen.ibis.ngs.vepsummary;

import java.util.HashSet;

public class Gene {
	
	private String id;
	private String symbol;
	private HashSet<String> transcript_ids;
	
	public Gene(String id) {
		this.id = id;
		this.symbol = "NA";
		this.transcript_ids = new HashSet<>();
	}
	
	public Gene(String id, String symbol) {
		this.id = id;
		this.symbol = symbol;
		this.transcript_ids = new HashSet<>();
	}
	
	public void addTranscript(String t) {
		this.transcript_ids.add(t);
	}
	
	public HashSet<String> getTranscripts() {
		return this.transcript_ids;
	}
	
	public String getSymbol()  {
		return this.symbol;
	}
	
	public String getId() {
		return this.id;
	}
	
	public String toString() {
		String res = id+"\t"+symbol+"\t";
		String transcripts = "";
		for(String t: transcript_ids) {
			transcripts += ","+t;
		}
		transcripts = transcripts.replaceFirst(",", "");
		
		return res+transcripts;
	}

}
package de.helmholtz_muenchen.ibis.ngs.lofsummary;

public class Annotation {
	
	private String gene_id;
	private String transcript_id;
	private String consequence;
	
	public Annotation(String gene_id, String transcript_id, String consequence) {
		this.gene_id = gene_id;
		this.transcript_id = transcript_id;
		this.consequence = consequence;
	}
	
	public String toString() {
		return "("+gene_id+"|"+transcript_id+"|"+consequence+")";
	}

}

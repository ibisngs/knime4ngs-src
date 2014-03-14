package de.helmholtz_muenchen.ibis.ngs.filtersummary.genes;

public class Gene {
	
	
	private String ensemblid;
	private String name;
	private String type;
	private int transcripts;
	private int afftrans;
	
	public Gene(String vat){
		
		String [] split = vat.split(":");
		
		if(split.length<7){
			System.out.println("Incorrect format: "+vat);
			return;
		}
		name=split[1];
		ensemblid = split[2];
		type = split[4];
		
		String [] trans = split[5].split("/");
		
		if(trans.length!=2){
			System.out.println("Incorrect transcript format: "+split[5]);
			return;
		}
		
		afftrans=Integer.parseInt(trans[0]);
		transcripts=Integer.parseInt(trans[1]);
		
	}
	
	public String toString(){
		return ensemblid+"\t"+name+"\t"+type+"\t"+transcripts+"\t"+afftrans;
	}
	
	public String getID(){
		return this.ensemblid;
	}
	
	public String getIDnoExt(){
		return this.ensemblid.split("\\.")[0];
	}
	
	public String getName(){
		return this.name;
	}
	
	public String getType(){
		return this.type;
	}
	
	public int getTranscripts(){
		return this.transcripts;
	}
	
	public int getAffTranscripts(){
		return this.afftrans;
	}
	
	

}

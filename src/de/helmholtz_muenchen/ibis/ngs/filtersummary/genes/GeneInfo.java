package de.helmholtz_muenchen.ibis.ngs.filtersummary.genes;

import java.util.HashMap;

import de.helmholtz_muenchen.ibis.utils.lofs.FileInputReader;

public class GeneInfo {
	
	private String geneinfofile;
	protected HashMap<String, String> ensembltomim;
	
	public GeneInfo(String gif){
		this.geneinfofile=gif;
		
		ensembltomim= new HashMap<String, String >(100000);
		
		loadGeneInfo();
	}
	
	private void loadGeneInfo(){
		
		FileInputReader f;
		try {
			f = new FileInputReader(geneinfofile);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Cannot load file: "+geneinfofile);
			return;
		}
		
		String line;
		
		while((line=f.read())!=null){
			
			if(line.startsWith("#")){
				continue;
			}
			
			// retrieve column with dbrefx
			String [] split=line.split("\t");
			
			// split db references
			String [] splitids=split[5].split("\\|");
			boolean ensembl=false;
			String eID="";
			boolean omim=false;
			String oID="";
			for(int i=0; i<splitids.length; i++){
				if(splitids[i].startsWith("Ensembl")){
					ensembl=true;
					eID=splitids[i].split(":")[1];
				}
				if(splitids[i].startsWith("MIM")){
					omim=true;
					oID=splitids[i].split(":")[1];
				}
			}
			//System.out.println(eID);
			//System.out.println(oID);
			
			//add id pair to hashmap
			if(omim && ensembl){
				if(ensembltomim.containsKey(eID)){
					//System.out.println(oID+"\t"+ensembltomim.get(eID));
					ensembltomim.put(eID, ensembltomim.get(eID)+","+oID);
				}
				else{
					ensembltomim.put(eID, oID);
				}
			}

		}
		
		try {
			f.closer();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public String getOmimID (String gene){		
		if(ensembltomim.containsKey(gene)){
			return ensembltomim.get(gene);
		}
		return "";
	}
	
//	public static void main (String args []){
//		GeneInfo gi= new GeneInfo("/home/marie-sophie/Uni/resources/Homo_sapiens.gene_info");
//	}

}

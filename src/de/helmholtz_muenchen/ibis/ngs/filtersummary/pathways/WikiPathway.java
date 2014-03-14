package de.helmholtz_muenchen.ibis.ngs.filtersummary.pathways;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;

import de.helmholtz_muenchen.ibis.utils.lofs.FileInputReader;

public class WikiPathway {
	
	private String pathwayfile;
	protected HashMap<String, LinkedList<String>> genetopathway;
	
	public WikiPathway(String file){
		this.pathwayfile=file;
		genetopathway= new HashMap<String, LinkedList<String>>(1000);
		
		loadPathways();
	}
	
	public HashMap<String, LinkedList<String>> getHashMap(){
		return this.genetopathway;
	}
	
	public String getPathway(String geneid){
		
		String res="";
		
		if(genetopathway.containsKey(geneid)){
			Iterator<String> j =genetopathway.get(geneid).iterator();
			boolean start=true;
			while(j.hasNext()){
				String geneID=j.next();
				if(start){
					start=false;
					res+=geneID;
				}
				else{
					res+=","+geneID;
				}
			}
		}
		
		return res;
	}
	
	private void loadPathways(){
		
		FileInputReader fir;
		try {
			fir = new FileInputReader(pathwayfile);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Cannot load file: "+genetopathway);
			return;
		}
		
		//skip first line
		String line=fir.read();
		
		while((line=fir.read())!=null){
			
			String [] split = line.split("\t");
			String pathway = split[0];
			
			if(split.length<10){
				continue;
			}
			String [] genes = split[9].split(",");
			for(int i =0; i<genes.length; i++){
				if(genetopathway.containsKey(genes[i])){
					genetopathway.get(genes[i]).add(pathway);
				}
				else{
					LinkedList<String> pathwaylist = new LinkedList<String>();
					pathwaylist.add(pathway);
					genetopathway.put(genes[i], pathwaylist);
				}
			}

		}
	}
	
	public void printHash(){
		Iterator<Entry<String, LinkedList<String>>> i = genetopathway.entrySet().iterator();
		
		while(i.hasNext()){
			Entry<String, LinkedList<String>> e = i.next();
			String res=e.getKey()+"\t";
			
			Iterator<String> j =e.getValue().iterator();
			boolean start=true;
			while(j.hasNext()){
				String geneID=j.next();
				if(start){
					start=false;
					res+=geneID;
				}
				else{
					res+=","+geneID;
				}
			}
			System.out.println(res);
		}
	}
	
//	public static void main (String args []){
//		WikiPathway wp = new WikiPathway("/home/marie-sophie/Uni/resources/wikipathways_data_Homo sapiens.tab");
//	}
}

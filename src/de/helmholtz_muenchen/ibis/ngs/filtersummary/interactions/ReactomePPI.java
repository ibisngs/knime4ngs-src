package de.helmholtz_muenchen.ibis.ngs.filtersummary.interactions;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;

import de.helmholtz_muenchen.ibis.utils.lofs.FileInputReader;

public class ReactomePPI {
	
	private String ppifile;
	private HashMap<String, LinkedList<String>> complexes;
	private HashMap<String, LinkedList<String>> reactions;
	private HashMap<String, LinkedList<String>> all;
	
	
	public ReactomePPI(String file){
		this.ppifile=file;
		complexes = new HashMap<String, LinkedList<String>>(200000);
		reactions = new HashMap<String, LinkedList<String>>(200000);
		all = new HashMap<String, LinkedList<String>>(200000);
		
		
		loadppis();
	}
	
	private void loadppis(){
		
		FileInputReader f;
		try {
			f = new FileInputReader(ppifile);
		} catch (Exception e) {
			System.out.println("Cannot read file: "+ppifile);
			return;
		}
		
		//skip first line
		String line=f.read();
		
		while((line=f.read())!=null){
			String[] splittab=line.split("\t");
			
			// no ensembl id available for one of the interaction partners
			if(splittab[1].equals("")|| splittab[4].equals("")){
				continue;
			}

			
			String [] eid1 = splittab[1].split("\\|");
			for(int x =0; x<eid1.length; x++){
				eid1[x]=eid1[x].split(":")[1];
			}
			String [] eid2 = splittab[4].split("\\|");
			for(int y=0; y<eid2.length; y++){
				eid2[y]=eid2[y].split(":")[1];
				
			}
			
			HashMap<String, LinkedList<String>> hashmap=null;
			// consider direct complex ppis
			if(splittab[6].equals("direct_complex")){
				hashmap=complexes;
			}
			if(splittab[6].equals("reaction")){
				hashmap=reactions;
			}
			
			if(hashmap!=null){
				
				for(int i = 0; i<eid1.length; i++){
					for(int j= 0; j<eid2.length; j++){
						//System.out.println(eid1[i]+"\t"+eid2[j]);
						
						if(hashmap.containsKey(eid1[i])){
							if(!hashmap.get(eid1[i]).contains(eid2[j])){
								hashmap.get(eid1[i]).add(eid2[j]);
							}
						}
						else{
							LinkedList<String> l= new LinkedList<String>();
							l.add(eid2[j]);
							hashmap.put(eid1[i], l);
						}
						
						if(hashmap.containsKey(eid2[j])){
							if(!hashmap.get(eid2[j]).contains(eid1[i])){
								hashmap.get(eid2[j]).add(eid1[i]);
							}
						}
						else{
							LinkedList<String> l= new LinkedList<String>();
							l.add(eid1[i]);
							hashmap.put(eid2[j], l);
							
						}
					}
				}
			}
			
			for(int i = 0; i<eid1.length; i++){
				for(int j= 0; j<eid2.length; j++){
					//System.out.println(eid1[i]+"\t"+eid2[j]);
					
					if(all.containsKey(eid1[i])){
						if(!all.get(eid1[i]).contains(eid2[j])){
							all.get(eid1[i]).add(eid2[j]);
						}
					}
					else{
						LinkedList<String> l= new LinkedList<String>();
						l.add(eid2[j]);
						all.put(eid1[i], l);
					}
					
					if(all.containsKey(eid2[j])){
						if(!all.get(eid2[j]).contains(eid1[i])){
							all.get(eid2[j]).add(eid1[i]);
						}
					}
					else{
						LinkedList<String> l= new LinkedList<String>();
						l.add(eid1[i]);
						all.put(eid2[j], l);
						
					}
				}
			}
		}		
	}
	
	public HashMap<String, LinkedList<String>> getHashMapComplex(){
		return this.complexes;
	}
	
	public HashMap<String, LinkedList<String>> getHashMapReaction(){
		return this.reactions;
	}
	
	public HashMap<String, LinkedList<String>> getHashMapAll(){
		return this.all;
	}
	
	public void printComplexes (){
		Iterator<Entry<String, LinkedList<String>>> it =complexes.entrySet().iterator();
		
		while(it.hasNext()){
			Entry<String, LinkedList<String>> e = it.next();
			System.out.println(e.getKey()+"\t"+e.getValue().size());
		}
		
	}

	
	public static void main (String [] args){
		ReactomePPI rp = new ReactomePPI("/home/marie-sophie/Uni/resources/homo_sapiens.interactions.txt");
		rp.printComplexes();
	}

}

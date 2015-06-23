package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class GenomeMap {
	
	private HashMap<String, ArrayList<ArrayList<Integer>>> genomeMap = new HashMap<>(24);

	/**
	 * @param args
	 */
	private String mapFile;
	
	public GenomeMap(String mapFile) {
		setMapFile(mapFile);
		createGenomeMap();
	}
	
	public String getMapFile() {
		return mapFile;
	}

	public void setMapFile(String mapFile) {
		this.mapFile = mapFile;
	}
		
	public HashMap<String, ArrayList<ArrayList<Integer>>> getGenomeMap() {
		return this.genomeMap;
	}
	
	
	private void createGenomeMap () {
		for (int i = 1; i < 23; i++) 
			this.genomeMap.put("chr"+i, new ArrayList<ArrayList<Integer>>());
		this.genomeMap.put("chrX", new ArrayList<ArrayList<Integer>>());
		this.genomeMap.put("chrY", new ArrayList<ArrayList<Integer>>());
		
		File input = new File (this.mapFile);
		
		try {
			Scanner sc = new Scanner(input, "UTF-8");
			
	        while (sc.hasNextLine()) {
	        	String currentLine = sc.nextLine().trim();
	        	
	        	String[] cols = currentLine.split("\t");
	        	/**
	        	 * Parsing the start and stop range for an N region
	        	 */
	        	int start = Integer.parseInt(cols[1]), stop = Integer.parseInt(cols[2]);
	    		ArrayList<Integer> range = new ArrayList<>(2);
	    		range.add(start);
	    		range.add(stop);
	    		/**
	    		 * parsing the chromosome for key
	    		 */
	        	String chr = cols[0];
	        	if (this.genomeMap.containsKey(chr)) {
	        		this.genomeMap.get(chr).add(range);
	        	}
	        	
	        }
//	        System.out.println(this.nMap.size());
//	        for (String s : this.nMap.keySet()) {
//	        	System.out.println(s + "\t" + this.nMap.get(s).size());
//	        }
	        
	    } 
	    catch (FileNotFoundException e) {
	        e.printStackTrace();
	    }
		
	}
	/*
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		GenomeMap n = new GenomeMap();
//		n.createNMap();
		
		for (String s : n.genomeMap.keySet())
			for (int i = 0; i < n.genomeMap.get(s).size(); i++)
			System.out.println(s + "\t" + n.genomeMap.get(s).get(i).get(0) + "\t" + n.genomeMap.get(s).get(i).get(1));
	}
*/
}

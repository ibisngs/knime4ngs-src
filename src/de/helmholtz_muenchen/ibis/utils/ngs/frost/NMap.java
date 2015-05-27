package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class NMap {

	/**
	 * @param args
	 */
	private final String file ="/home/ibis/tanzeem.haque/Documents/hg19/hg19_N_regions_tab.txt";
	private HashMap<String, ArrayList<ArrayList<Integer>>> nMap = new HashMap<>(24);
	
	public HashMap<String, ArrayList<ArrayList<Integer>>> getNMap() {
	
		createNMap();
//		for (String s : this.nMap.keySet())
//			for (int i = 0; i < this.nMap.get(s).size(); i++)
//				System.out.println(s + "\t" + this.nMap.get(s).get(i).get(0) + "\t" + this.nMap.get(s).get(i).get(1));

		return this.nMap;
	}
	
	
	private void createNMap () {
		for (int i = 1; i < 23; i++) 
			this.nMap.put("chr"+i, new ArrayList<ArrayList<Integer>>());
		this.nMap.put("chrX", new ArrayList<ArrayList<Integer>>());
		this.nMap.put("chrY", new ArrayList<ArrayList<Integer>>());
		
		File input = new File (this.file);
		
		try {
			Scanner sc = new Scanner(input, "UTF-8");
			
	        while (sc.hasNextLine()) {
	        	String currentLine = sc.nextLine().trim();
//	        	if(currentLine == null )
//		            throw new IOException( this.file + " File is empty " );  
	        	
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
	        	if (this.nMap.containsKey(chr)) {
//	        		System.out.println(chr + "\t" + this.nMap.get(chr).size());
	        		this.nMap.get(chr).add(range);
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
		NMap n = new NMap();
		n.createNMap();
		
//		for (String s : n.nMap.keySet())
//			for (int i = 0; i < n.nMap.get(s).size(); i++)
//			System.out.println(s + "\t" + n.nMap.get(s).get(i).get(0) + "\t" + n.nMap.get(s).get(i).get(1));
	}
*/
}

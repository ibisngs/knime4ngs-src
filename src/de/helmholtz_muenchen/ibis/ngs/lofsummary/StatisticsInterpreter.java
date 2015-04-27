package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.TreeMap;

public class StatisticsInterpreter {
	
	public static void main (String [] args) throws IOException {
		
		String file = "/home/ibis/tim.jeske/KORAdata/analysis_ready/20150413_loftee_gmaf_sift_polyphen/hc_lof_stats/variant_effect_LoF.sample_statistic.txt";
		HashSet<String> myGenes = new HashSet<>();
		ArrayList<String> sample_ids = new ArrayList<>();
		TreeMap<String, Integer> gene2occ = new TreeMap<>();
		
		//sample id -> gene id -> true/false
		HashMap<String, HashMap<String, String>> table = new HashMap<>();
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		br.readLine(); //ignore header
		String line;
		int samples = 0;
		while((line=br.readLine())!=null) {
			samples++;
			String [] fields = line.split("\t");
			String [] genes = fields[5].split(",");
			String sample_id = fields[0];
			sample_ids.add(sample_id);
			for(String s:genes) {
				String gene_id = s.split(":")[0];
				String gene_symbol = s.split(":")[1];
//				if(gene_symbol.startsWith("OR")) continue;
//				if(gene_symbol.startsWith("KRTA")) continue;
				myGenes.add(gene_id);

				HashMap<String, String> tmp = table.get(sample_id);
				if(tmp==null) {
					tmp = new HashMap<>();
				}
				
				tmp.put(gene_id, "1");
				table.put(sample_id, tmp);

				Integer i = gene2occ.get(gene_id);
				if(i==null) i = 0;
				i++;
				gene2occ.put(gene_id, i);
				
			}
		}
		br.close();
		
		for(String s: myGenes) {
			System.out.println(s);
		}
		
		System.out.println(myGenes.size()+" genes are completely disabled.");
		System.out.println("#############");
		
		for(Entry<String, Integer> entry:gene2occ.entrySet()) {
			double percentage = ((double)entry.getValue())/samples;
			System.out.println(entry.getKey() +"\t"+percentage);
		}
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(file.replace("txt","table.csv")));
		bw.write("sample_id");
		for(String g:myGenes) {
			bw.write("\t"+g);
		}
		bw.newLine();
		
		for(String s: sample_ids) {
			bw.write(s);
			for(String g:myGenes ) {
				String val = table.get(s).get(g);
				if(val == null) val = "0";
				bw.write("\t"+val);
			}
			bw.newLine();
		}
		bw.close();
		
	}

}

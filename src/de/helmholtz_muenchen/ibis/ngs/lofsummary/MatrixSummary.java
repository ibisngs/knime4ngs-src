package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;

import de.helmholtz_muenchen.ibis.utils.ngs.ContingencyTable;

public class MatrixSummary {

	private String [] samples;
	private HashMap<String, HashMap<Integer,Short>> transcript2sample2aff;
	
	public MatrixSummary(String file) throws IOException{
		
		this.transcript2sample2aff = new HashMap<>();
		
		BufferedReader br = Files.newBufferedReader(Paths.get(file));
		String [] header = br.readLine().split("\t");
		
		samples = new String[header.length-1];
		for(int i = 1; i < header.length; i++) {
			samples[i-1] = header[i];
		}
		
		String line,t;
		String [] fields;
		HashMap<Integer,Short> tmp;
		while((line=br.readLine())!=null) {
			tmp = new HashMap<>();
			fields = line.split("\t");
			t = fields[0];
			for(int i = 1; i < fields.length; i++) {
				if(fields[i].equals("unaff")) continue;
				tmp.put(i-1, replaceAff(fields[i]));
			}
			transcript2sample2aff.put(t, tmp);
		}
		br.close();
	}
	
	private Short replaceAff(String aff) {
		if(aff.equals("het") || aff.equals("pot_comp")) {
			return 1;
		} else if(aff.equals("hom")) {
			return 2;
		} 
		return -1;
	}
	
	/**
	 * Converts matrix to a contingency table for each gene
	 * @return a map with genes as keys (id_symbol) and ContingencyTables as values
	 */
	public HashMap<String, ContingencyTable> toTables() {
		HashMap<String, ContingencyTable> result = new HashMap<>();
		
		HashMap<String, HashSet<Integer>> gene2aff_cases = new HashMap<>();
		HashMap<String, HashSet<Integer>> gene2un_cases = new HashMap<>();
		HashMap<String, HashSet<Integer>> gene2aff_ctrls = new HashMap<>();
		HashMap<String, HashSet<Integer>> gene2un_ctrls = new HashMap<>();
		HashSet<Integer> cases = new HashSet<>();
		HashSet<Integer> ctrls = new HashSet<>();
		
		for(int i = 0; i< samples.length;i++) {
			if(samples[i].contains("case")) {
				cases.add(i);
			} else if(samples[i].contains("control")) {
				ctrls.add(i);
			}
		}
		
		String gene;
		Short aff;
		HashSet<Integer> tmp;
		HashMap<Integer,Short> sample2aff;
		for(String transcript : transcript2sample2aff.keySet()) {
			gene = transcript.split("_",2)[1]; //gene = id_symbol
			if(!gene2aff_ctrls.containsKey(gene)) {
				gene2aff_ctrls.put(gene, new HashSet<>());
			}
			if(!gene2aff_cases.containsKey(gene)) {
				gene2aff_cases.put(gene, new HashSet<>());
			}
			if(!gene2un_ctrls.containsKey(gene)) {
				gene2un_ctrls.put(gene, new HashSet<>(ctrls));
			}
			if(!gene2un_cases.containsKey(gene)) {
				gene2un_cases.put(gene, new HashSet<>(cases));
			}
			
			sample2aff = transcript2sample2aff.get(transcript);
			
			for(Integer s: sample2aff.keySet()) {
				aff = sample2aff.get(s);
				if(aff >= 1) { //affected
					if(samples[s].contains("case")) {
						tmp = gene2aff_cases.get(gene);
						tmp.add(s);
						gene2aff_cases.put(gene, tmp);
						gene2un_cases.get(gene).remove(s);
					} else if(samples[s].contains("control")) {
						tmp = gene2aff_ctrls.get(gene);
						tmp.add(s);
						gene2aff_ctrls.put(gene, tmp);
						gene2un_ctrls.get(gene).remove(s);
					}
				} else { //unclear whether affected or not
					if(samples[s].contains("case")) {
						gene2un_cases.get(gene).remove(s);
					} else if(samples[s].contains("control")) {
						gene2un_ctrls.get(gene).remove(s);
					}
				}
			}
		}
		
		int case_aff, control_aff;
		for(String g: gene2un_ctrls.keySet()) {
			case_aff = gene2aff_cases.get(g).size();
			control_aff = gene2aff_ctrls.get(g).size();
			if(case_aff == 0 && control_aff == 0) {
    			continue;
    		}
			result.put(g, new ContingencyTable(case_aff,gene2un_cases.get(g).size(),control_aff,gene2un_ctrls.get(g).size()));
		}
		return result;
	}
	
	
}

package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

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
	public HashMap<String, ContingencyTable> toTables(Identifier identifier) {
		HashMap<String, ContingencyTable> result = new HashMap<>();
		
		HashMap<String, HashSet<Integer>> entity2aff_cases = new HashMap<>();
		HashMap<String, HashSet<Integer>> entity2un_cases = new HashMap<>();
		HashMap<String, HashSet<Integer>> entity2aff_ctrls = new HashMap<>();
		HashMap<String, HashSet<Integer>> entity2un_ctrls = new HashMap<>();
		HashSet<Integer> cases = new HashSet<>();
		HashSet<Integer> ctrls = new HashSet<>();
		
		for(int i = 0; i< samples.length;i++) {
			if(samples[i].contains("case")) {
				cases.add(i);
			} else if(samples[i].contains("control")) {
				ctrls.add(i);
			}
		}
		
		List<String> entities;
		Short aff;
		HashMap<Integer,Short> sample2aff;
		for(String transcript : transcript2sample2aff.keySet()) {
			entities = identifier.getMappings(transcript);
//			gene = transcript.split("_",2)[1]; //gene = id_symbol
			for(String entity: entities) {
				if(!entity2aff_cases.containsKey(entity)) {
					entity2aff_cases.put(entity, new HashSet<>());
					entity2aff_ctrls.put(entity, new HashSet<>());
					entity2un_ctrls.put(entity, new HashSet<>(ctrls));
					entity2un_cases.put(entity, new HashSet<>(cases));
				}
				
				sample2aff = transcript2sample2aff.get(transcript);
				
				for(Integer s: sample2aff.keySet()) {
					aff = sample2aff.get(s);
					if(aff >= 1) { //affected
						if(samples[s].contains("case")) {
							entity2aff_cases.get(entity).add(s);
							entity2un_cases.get(entity).remove(s);
						} else if(samples[s].contains("control")) {
							entity2aff_ctrls.get(entity).add(s);
							entity2un_ctrls.get(entity).remove(s);
						}
					} else { //unclear whether affected or not
						if(samples[s].contains("case")) {
							entity2un_cases.get(entity).remove(s);
						} else if(samples[s].contains("control")) {
							entity2un_ctrls.get(entity).remove(s);
						}
					}
				}
			}
		}
		
		int case_aff, control_aff;
		for(String e: entity2un_ctrls.keySet()) {
			case_aff = entity2aff_cases.get(e).size();
			control_aff = entity2aff_ctrls.get(e).size();
			if(case_aff == 0 && control_aff == 0) {
    			continue;
    		}
			result.put(e, new ContingencyTable(case_aff,entity2un_cases.get(e).size(),control_aff,entity2un_ctrls.get(e).size()));
		}
		return result;
	}
}

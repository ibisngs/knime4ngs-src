package de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import de.helmholtz_muenchen.ibis.ngs.geneBasedAnalysis.CaseControlArray;
import de.helmholtz_muenchen.ibis.utils.ngs.ContingencyTable;

public class MatrixSummary {

	private String [] samples;
	private HashMap<String, HashMap<Integer,Short>> transcript2sample2aff;
	HashSet<Integer> cases,ctrls; //contain indicies of samples array
	
	public MatrixSummary(String file) throws IOException{
		
		this.transcript2sample2aff = new HashMap<>();
		
		BufferedReader br = Files.newBufferedReader(Paths.get(file));
		String [] header = br.readLine().split("\t");
		
		samples = new String[header.length-1];
		for(int i = 1; i < header.length; i++) {
			samples[i-1] = header[i];
		}
		
		cases = new HashSet<>();
		ctrls = new HashSet<>();
		for(int i = 0; i< samples.length;i++) {
			if(samples[i].contains("case")) {
				cases.add(i);
			} else if(samples[i].contains("control")) {
				ctrls.add(i);
			}
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
		if(aff.equals("het") || aff.equals("pot_comp") || aff.equals("mat") || aff.equals("pat")) {
			return 1;
		} else if(aff.equals("hom") || aff.equals("comp")) {
			return 2;
		} 
		return -1;
	}
	
	/**
	 * Converts matrix to a contingency table for each gene
	 * @return a map with keys defined by the Identifier implementation and ContingencyTables as values
	 */
	public HashMap<String, ContingencyTable> toTables(Identifier identifier) {
		HashMap<String, ContingencyTable> result = new HashMap<>();
		
		HashMap<String, HashSet<Integer>> entity2aff_cases = new HashMap<>();
		HashMap<String, HashSet<Integer>> entity2un_cases = new HashMap<>();
		HashMap<String, HashSet<Integer>> entity2aff_ctrls = new HashMap<>();
		HashMap<String, HashSet<Integer>> entity2un_ctrls = new HashMap<>();
		
		List<String> entities;
		Short aff;
		HashMap<Integer,Short> sample2aff;
		for(String transcript : transcript2sample2aff.keySet()) {
			sample2aff = transcript2sample2aff.get(transcript);
			entities = identifier.getMappings(transcript);
			for(String entity: entities) {
				if(!entity2aff_cases.containsKey(entity)) {
					entity2aff_cases.put(entity, new HashSet<>());
					entity2aff_ctrls.put(entity, new HashSet<>());
					entity2un_ctrls.put(entity, new HashSet<>(ctrls));
					entity2un_cases.put(entity, new HashSet<>(cases));
				}
				
				for(Integer s: sample2aff.keySet()) {
					aff = sample2aff.get(s);
					if(samples[s].contains("case")) {
						if(aff >= 1) {
							entity2aff_cases.get(entity).add(s);
						}
						entity2un_cases.get(entity).remove(s);
					} else if(samples[s].contains("control")) {
						if(aff >=1) {
							entity2aff_ctrls.get(entity).add(s);
						}
						entity2un_ctrls.get(entity).remove(s);
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
	
	/**
	 * simplified implementation of toTables(Identifier identifier) used for analyses on transcript level
	 * @return
	 */
	public HashMap<String, ContingencyTable> toTables() {
		HashMap<String, ContingencyTable> result = new HashMap<>();
		
		int aff_cases, un_cases, aff_ctrls, un_ctrls;
		
		Short aff;
		HashMap<Integer,Short> sample2aff;
		for(String transcript : transcript2sample2aff.keySet()) {
			
			aff_cases = 0;
			aff_ctrls = 0;
			un_cases = cases.size();
			un_ctrls = ctrls.size();
			sample2aff = transcript2sample2aff.get(transcript);
				
			for(Integer s: sample2aff.keySet()) {
				aff = sample2aff.get(s);
				 //affected
				if(samples[s].contains("case")) {
					if(aff >= 1) {
						aff_cases++;
					}
					un_cases--;
				} else if(samples[s].contains("control")) {
					if(aff >=1) {
						aff_ctrls++;
					}		
					un_ctrls--;
				}
			}
			if(aff_cases == 0 && aff_ctrls == 0) {
				continue;
			}
			result.put(transcript, new ContingencyTable(aff_cases,un_cases,aff_ctrls,un_ctrls));
		}
		return result;
	}
	
	public HashMap<String,CaseControlArray> toArrays() {
		HashMap<String,CaseControlArray> res = new HashMap<>();
		
		short [] cases,controls;
		Short aff;
		HashMap<Integer,Short> sample2aff;
		ArrayList<Short> my_cases, my_ctrls;
		int un_case, un_ctrl;
		
		for(String transcript : transcript2sample2aff.keySet()) {
			
			un_case = this.cases.size();
			my_cases = new ArrayList<Short>();
			un_ctrl = this.ctrls.size();
			my_ctrls = new ArrayList<Short>();
			
			sample2aff = transcript2sample2aff.get(transcript);
				
			for(Integer s: sample2aff.keySet()) {
				aff = sample2aff.get(s);
				 //affected
				if(samples[s].contains("case")) {
					if(aff >= 1) {
						my_cases.add(aff);
					}
					un_case--;
				} else if(samples[s].contains("control")) {
					if(aff>= 1) {
						my_ctrls.add(aff);
					}
					un_ctrl--;
				}
			}
			if(my_cases.size() == 0 && my_ctrls.size() == 0) {
				continue;
			}
			cases = new short[my_cases.size()+un_case];
			for(int i = 0; i < my_cases.size(); i++) {
				cases[i] = my_cases.get(i);
			}
			controls = new short[my_ctrls.size()+un_ctrl];
			for(int i = 0; i < my_ctrls.size(); i++) {
				controls[i] = my_ctrls.get(i);
			}
			res.put(transcript, new CaseControlArray(cases,controls));
		}
		
		return res;
	}
	
	public HashMap<String,CaseControlArray> toArrays(Identifier identifier) {
		HashMap<String,CaseControlArray> res = new HashMap<>();
		
		HashMap<String, HashSet<Integer>> entity2hom_cases = new HashMap<>();
		HashMap<String, HashSet<Integer>> entity2het_cases = new HashMap<>();
		HashMap<String, HashSet<Integer>> entity2un_cases = new HashMap<>();
		HashMap<String, HashSet<Integer>> entity2hom_ctrls = new HashMap<>();
		HashMap<String, HashSet<Integer>> entity2het_ctrls = new HashMap<>();
		HashMap<String, HashSet<Integer>> entity2un_ctrls = new HashMap<>();
		
		List<String> entities;
		Short aff;
		HashMap<Integer,Short> sample2aff;
		
		for(String transcript : transcript2sample2aff.keySet()) {
			sample2aff = transcript2sample2aff.get(transcript);
			entities = identifier.getMappings(transcript);
			for(String entity: entities) {
				if(!entity2hom_cases.containsKey(entity)) {
					entity2hom_cases.put(entity, new HashSet<>());
					entity2hom_ctrls.put(entity, new HashSet<>());
					entity2het_cases.put(entity, new HashSet<>());
					entity2het_ctrls.put(entity, new HashSet<>());
					entity2un_ctrls.put(entity, new HashSet<>(ctrls));
					entity2un_cases.put(entity, new HashSet<>(cases));
				}
				
				for(Integer s: sample2aff.keySet()) {
					aff = sample2aff.get(s);
					if(samples[s].contains("case")) {
						if(aff==1) {
							entity2het_cases.get(entity).add(s);
						} else if (aff==2) {
							entity2hom_cases.get(entity).add(s);
						}
						entity2un_cases.get(entity).remove(s);
					} else if(samples[s].contains("control")) {
						if(aff==1) {
							entity2het_ctrls.get(entity).add(s);
						} else if (aff==2) {
							entity2hom_ctrls.get(entity).add(s);
						}
						entity2un_ctrls.get(entity).remove(s);
					}
				}
			}
		}

		ArrayList<Short> my_cases,my_controls;
		for(String e: entity2un_ctrls.keySet()) {
			if(entity2hom_cases.size() == 0 && entity2het_cases.size() == 0 && entity2hom_ctrls.size() == 0 && entity2het_ctrls.size() == 0) {
    			continue;
    		}
			my_cases = new ArrayList<>();
			my_controls = new ArrayList<>();
			for(int i = 0; i < entity2hom_cases.get(e).size(); i++) {
				my_cases.add((short)2);
			}
			for(int i = 0; i < entity2het_cases.get(e).size(); i++) {
				my_cases.add((short)1);
			}
			for(int i = 0; i < entity2un_cases.get(e).size(); i++) {
				my_cases.add((short)0);
			}
			for(int i = 0; i < entity2hom_ctrls.get(e).size(); i++) {
				my_controls.add((short)2);
			}
			for(int i = 0; i < entity2het_ctrls.get(e).size(); i++) {
				my_controls.add((short)1);
			}
			for(int i = 0; i < entity2un_ctrls.get(e).size(); i++) {
				my_controls.add((short)0);
			}
			res.put(e, new CaseControlArray(my_cases,my_controls));
		
		}
		return res;
	}
	
	public HashMap<String,Double> getFrequencies (Identifier identifier) {
		
		HashMap<String, Double> result = new HashMap<>();
		
		HashMap<String, HashSet<Integer>> entity2aff = new HashMap<>();
		HashMap<String, HashSet<Integer>> entity2un = new HashMap<>();
		
		List<String> entities;
		Short aff;
		HashSet<Integer> my_samples = new HashSet<>();
		for(int i = 0; i < samples.length; i++) {
			my_samples.add(i);
		}
		
		HashMap<Integer,Short> sample2aff;
		for(String transcript : transcript2sample2aff.keySet()) {
			sample2aff = transcript2sample2aff.get(transcript);
			entities = identifier.getMappings(transcript);
			
			for(String entity: entities) {
				if(!entity2aff.containsKey(entity)) {
					entity2aff.put(entity, new HashSet<>());
					entity2un.put(entity, new HashSet<>(my_samples));
				}
				
				for(Integer s: sample2aff.keySet()) {
					aff = sample2aff.get(s);
					if(aff >= 1) { //affected
						entity2aff.get(entity).add(s);
					}
					entity2un.get(entity).remove(s);
				}
			}
		}
		
		int affected, unaff;
		for(String e: entity2un.keySet()) {
			affected = entity2aff.get(e).size();
			unaff = entity2un.get(e).size();
			result.put(e, (double)affected/(double)(affected+unaff));
		}
		return result;
	}
	
	public HashMap<String,Double> getFrequencies () {
		
		HashMap<String, Double> result = new HashMap<>();
		
		int affected, unaff;
		Short aff;
		
		HashMap<Integer,Short> sample2aff;
		for(String transcript : transcript2sample2aff.keySet()) {
			affected = 0;
			unaff = samples.length;
				
			sample2aff = transcript2sample2aff.get(transcript);
				
			for(Integer s: sample2aff.keySet()) {
				aff = sample2aff.get(s);
				if(aff >= 1) { //affected
					affected++;
				}
				unaff--;
			}
			
			result.put(transcript.split("_")[0], (double)affected/(double)(affected+unaff));
		}
		return result;
	}
}
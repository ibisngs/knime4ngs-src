package de.helmholtz_muenchen.ibis.ngs.extendgenesummary;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class StringGraph {
	
	private static HashMap<String,Values> gene_values = new HashMap<>();
	private static HashMap<String,String> protein_gene = new HashMap<>();
	private static HashMap<String,HashSet<String>> protein_graph = new HashMap<>();
	
	
	public static void readGraph(String file) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = br.readLine(); //omit header
		String [] fields;
		while((line=br.readLine())!=null) {
			fields = line.split(" ");
			if(Integer.parseInt(fields[2])>400) {
				String pepA = fields[0].split("\\.")[1];
				String pepB = fields[1].split("\\.")[1];
				if (protein_graph.containsKey(pepA)) {
					protein_graph.get(pepA).add(pepB);
				} else {
					HashSet<String> tmp = new HashSet<>();
					tmp.add(pepB);
					protein_graph.put(pepA, tmp);
				}

				if (protein_graph.containsKey(pepB)) {
					protein_graph.get(pepB).add(pepA);
				} else {
					HashSet<String> tmp = new HashSet<>();
					tmp.add(pepA);
					protein_graph.put(pepB, tmp);
				}
			}
		}
	}
	
	public static void readPepFa(String file) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		String [] fields;
		String pep;
		String gene;
		while((line=br.readLine())!=null) {
			if(line.startsWith(">")) {
				fields = line.split(" ");
				pep = fields[0].replaceFirst(">", "");
				gene = fields[3].split(":")[1];
				protein_gene.put(pep, gene);
			}
			
		}
		br.close();
	}
	
	public static void readGeneSummary(String file) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = br.readLine();
		String [] fields;
		String gene;
		int lofcount;
		double p_lof;
		while((line=br.readLine())!=null) {
			fields = line.split("\t");
			gene = fields[0];
			lofcount = Integer.parseInt(fields[2]) + Integer.parseInt(fields[3]);
			p_lof = Double.parseDouble(fields[4]);
			gene_values.put(gene, new Values(lofcount,p_lof,0,0.0));
		}
		
	}
	
	public static void calcDegree() {
		for(String pep:protein_graph.keySet()) {
			if(protein_gene.containsKey(pep)) {
				String gene = protein_gene.get(pep);
				int degree = protein_graph.get(pep).size();
				if(gene_values.containsKey(gene)) {
					Values v = gene_values.get(gene);
					v.setMax_degree(degree);
				}
			}
		}
	}
	
	private static boolean isEdge(String a, String b) {
		if(protein_graph.containsKey(a) && protein_graph.containsKey(b)) {
			return protein_graph.get(a).contains(b);
		}
		return false;
	}
	
	public static void calcClusterCoeff() {
		for(String pep:protein_graph.keySet()) {
			if(protein_gene.containsKey(pep)) {
				String gene = protein_gene.get(pep);
				if(gene_values.containsKey(gene)) {
					Values v = gene_values.get(gene);
					
					Object [] neighbours = protein_graph.get(pep).toArray();
					int possible_edges = (int) (neighbours.length*(neighbours.length-1)*0.5);
					int existing_edges = 0;
					for(int i = 0; i< neighbours.length-1; i++) {
						if(isEdge((String)neighbours[i],(String)neighbours[i+1])) {
							existing_edges++;
						}
					}
					v.setMax_cc((double)existing_edges/(double)possible_edges);
				}
			}
		}
	}
	
	public static void writeResult(String file) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		bw.write("gene_id\tlofcount\tp_lof\tmax_deg\tmax_cc");
		bw.newLine();
		for(String gene: gene_values.keySet()) {
			Values v = gene_values.get(gene);
			bw.write(gene+"\t");
			bw.write(v.getLof_count()+"\t");
			bw.write(v.getP_lof()+"\t");
			bw.write(v.getMax_degree()+"\t");
			bw.write(v.getMax_cc()+"");
			bw.newLine();
		}
		
		bw.close();
		
	}
	
	public static void main (String [] args) throws IOException {
		String gene_sum = "/home/tim/Dropbox/Studium/Helmholtz/Master/ExAC_Analysis/filtered/hcLoFs/ExAC.r0.3.PASS.vep.filtered.gene_summary.tsv";
		String pep_fa = "/home/tim/LOF_Project_local/files/Homo_sapiens.GRCh37.75.pep.all.fa";
		String string_graph = "/home/tim/LOF_Project_local/files/9606.protein.links.v10.txt";
		
		System.out.println("Read gene summary");
		readGeneSummary(gene_sum);
		System.out.println("Read pep fasta");
		readPepFa(pep_fa);
		System.out.println("Read graph");
		readGraph(string_graph);
		
		System.out.println("Calc degree");
		calcDegree();
		System.out.println("Calc cluster coefficient");
		calcClusterCoeff();
		
		System.out.println("Write result");
		writeResult("/home/tim/Dropbox/Studium/Helmholtz/Master/ExAC_Analysis/filtered/hcLoFs/extended_gene_summary.tsv");
	}

}

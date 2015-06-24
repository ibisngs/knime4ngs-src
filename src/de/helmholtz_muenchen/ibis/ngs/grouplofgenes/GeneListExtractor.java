package de.helmholtz_muenchen.ibis.ngs.grouplofgenes;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.knime.core.node.NodeLogger;

public class GeneListExtractor {
	
	HashSet<String> genes;
	HashMap<String, Double> gene_prob;
	int samples;
	String gene_summary;
	
	final static NodeLogger logger = NodeLogger.getLogger(GroupLoFGenesNodeModel.class);
	
	public GeneListExtractor(String cds_file, String geneSummary, int samples) {
		try {
			this.genes = readCDSFile(cds_file);
		} catch (IOException e) {
			logger.warn(cds_file+" could not be read. Never observed genes cannot be found.");
		}
		try {
			this.gene_prob = readGeneSummary(geneSummary);
		} catch (IOException e) {
			logger.error(geneSummary +" could not be read.");
		}
		this.gene_summary = geneSummary;
		this.samples = samples;
	}
	
	public String [] getGeneLists() {
		String [] outfiles = new String[3];
		
		HashSet<String> never_observed = new HashSet<>();
		HashSet<String> lof_intolerant = new HashSet<>();
		HashSet<String> lof_tolerant = new HashSet<>();
		
		double avg_prob = 0.0;
		for(String s:gene_prob.keySet()) {
			avg_prob += gene_prob.get(s);
		}
		avg_prob = avg_prob/gene_prob.size();
		logger.info("average probability: "+avg_prob);
		
		Set<String> all_genes;
		if(genes.size()==0) {
			all_genes = gene_prob.keySet();
		} else {
			all_genes = genes;
		}
		
		double prob;
		int observed;
		BinomialDistribution bd= new BinomialDistribution(samples,avg_prob);
		double p_value;
		
		
		for(String s:all_genes) {
			if(gene_prob.containsKey(s)) {
				prob = gene_prob.get(s);
				observed = (int) (prob*samples);
				if(prob > avg_prob) {
					p_value = 1 - bd.cumulativeProbability(observed-1);
					if(p_value < 0.05) {
						lof_tolerant.add(s);
					}
				} else if (prob < avg_prob) {
					p_value = bd.cumulativeProbability(observed);
					if(p_value < 0.05) {
						lof_intolerant.add(s);
					}
				}
				
			} else {
				never_observed.add(s);
			}
		}
		
		String never_obs = gene_summary.replace("gene_summary", "never_observed_genes");
		try {
			this.writeGeneList(never_observed, never_obs);
		} catch (IOException e) {
			logger.error(never_obs+" could not be written: "+e.getMessage());
		}
		outfiles[0] = never_obs;
		
		String lof_tol = gene_summary.replace("gene_summary", "lof_tolerant_genes");
		try {
			this.writeGeneList(lof_tolerant, lof_tol);
		} catch (IOException e) {
			logger.error(lof_tol+" could not be written: "+e.getMessage());
		}
		outfiles[1] = lof_tol;
		
		String lof_in = gene_summary.replace("gene_summary", "lof_intolerant_genes");
		try {
			this.writeGeneList(lof_intolerant, lof_in);
		} catch (IOException e) {
			logger.error(lof_in+" could not be written: "+e.getMessage());
		}
		outfiles[2] = lof_in;
		
		return outfiles;
	}
	
	private HashSet<String> readCDSFile(String cds_file) throws IOException {
		HashSet<String> result = new HashSet<>();
		String [] fields;
		String gene_id;
		
		FileInputStream inputStream = null;
		Scanner sc = null;
		try {
		    inputStream = new FileInputStream(cds_file);
		    sc = new Scanner(inputStream, "UTF-8");
		    while (sc.hasNextLine()) {
		        String line = sc.nextLine();
		        if(line.startsWith(">")) {
		        	fields = line.split("\\s");
		        	gene_id = fields[3].split(":")[1];
		        	result.add(gene_id);
		        }
		    }
		    // note that Scanner suppresses exceptions
		    if (sc.ioException() != null) {
		        throw sc.ioException();
		    }
		} finally {
		    if (inputStream != null) {
		        inputStream.close();
		    }
		    if (sc != null) {
		        sc.close();
		    }
		}
		return result;
	}
	
	private HashMap<String,Double> readGeneSummary(String file) throws IOException {
		HashMap<String,Double> result = new HashMap<>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		//ignore headline
		String line = br.readLine();
		String[] fields;
		while((line = br.readLine())!=null) {
			fields = line.split("\t");
			result.put(fields[0], Double.parseDouble(fields[4]));
			
		}
		br.close();
		return result;
	}
	
	private void writeGeneList(HashSet<String> set,String file) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		for(String s:set) {
			bw.write(s);
			bw.newLine();
		}
		bw.close();
	}
}

package de.helmholtz_muenchen.ibis.ngs.annotationComparator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;


public class Comparison {

	private final static String NEWLINE = System.getProperty("line.separator");
	
	private ArrayList<StatLoF> first, second;
	
	/**
	 * tool A/ B	full	partial	no_variant
	 * full			...		...		...
	 * partial		...		...		...
	 * no_variant	...		...		...
	 */
	private int[][] effectTab;
	
	/** 
	 * tool A/B				stop_gained	splice_site_variant	frameshift_variant	no_variant
	 * stop_gained			...			...					...					...
	 * splice_site_variant	...			...					...					...
	 * frameshift_variant	...			...					...					...
	 * no_variant			...			...					...					...
	 */
	private int[][] consequenceTab;
	
	public Comparison(String file1, String file2) {
		try {
			first = readStatFile(file1);
			second = readStatFile(file2);
			effectTab = new int[3][3];
			consequenceTab = new int[4][4];
		} catch (IOException e) {
			System.err.println(file1+ " could not be read!");
			e.printStackTrace();
		}
	}
	
	private ArrayList<StatLoF> readStatFile(String file) throws IOException {
		
		String line, chr, pos, ref_allele, alt_allele;
		String [] fields, consequences, effects;
		List<String> header;
		
		ArrayList<StatLoF> result = new ArrayList<>();
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		header = Arrays.asList(br.readLine().split("\t")); //header
		while((line=br.readLine())!=null) {
			fields = line.split("\t");
			chr = fields[header.indexOf("chr")];
			pos = fields[header.indexOf("pos")];
			ref_allele = fields[header.indexOf("ref_allele")];
			alt_allele = fields[header.indexOf("alt_allele")];
			consequences = fields[header.indexOf("consequence")].split(",");
			effects = fields[header.indexOf("effect")].split(",");
//			System.out.println(chr+" "+pos+" "+alt_allele+" "+consequences[0]+" "+effects[0]);
			result.add(new StatLoF(chr, pos, ref_allele, alt_allele, consequences, effects));
		}
		return result;
	}
	
	public void calcDiff() {
		int i = 0;
		int j = 0;
		int comp = -1;
		while (i < first.size() && j < second.size()) {
			StatLoF lof1 = first.get(i);
			StatLoF lof2 = second.get(j);
			comp = lof1.comparePos(lof2);

			if (comp == -1) {
				i++;
			} else if (comp == 1) {
				j++;
			} else {
				//create two lists and compare each
				ArrayList<StatLoF> tmp1 = new ArrayList<>();
				tmp1.add(lof1);
				ArrayList<StatLoF> tmp2 = new ArrayList<>();
				tmp2.add(lof2);
				i++;
				j++;
				while(i< first.size() && first.get(i).comparePos(lof1)==0) {
					tmp1.add(first.get(i));
					i++;
				}
				while(j < second.size() && second.get(j).comparePos(lof2)==0) {
					tmp2.add(second.get(j));
					j++;
				}
				for(int k = 0; k < tmp1.size(); k++) {
					for(int l = k; l < tmp2.size(); l++) {
						StatLoF mylof1 = tmp1.get(k);
						StatLoF mylof2 = tmp2.get(l);
//						System.out.println(mylof1);
//						System.out.println(mylof2);
//						System.out.println("++++++++++++++++++++++++++");
						HashSet<String> tmp_cons1 = mylof1.getConsequences();
						HashSet<String> tmp_cons2 = mylof2.getConsequences();
						HashSet<String> tmp_eff1 = mylof1.getEffects();
						HashSet<String> tmp_eff2 = mylof2.getEffects();
						
						//fill effect table

					}
				}
			}
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String file1,file2;
		
		file1 = "/home/ibis/tim.jeske/KORAdata/analysis_ready/20150319_vat_k/analysis_ready.diabetes.filtered.haplotypecaller.vat.lof_statistic.txt";
		file2 = "/home/ibis/tim.jeske/KORAdata/analysis_ready/20150401_loftee_q/all_lof_stats/variant_effect_LoF.lof_statistic.txt";
		
		Comparison my = new Comparison(file1,file2);
		my.calcDiff();

	}

}

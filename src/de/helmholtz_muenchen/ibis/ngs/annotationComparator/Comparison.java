package de.helmholtz_muenchen.ibis.ngs.annotationComparator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class Comparison {

	private final static String NEWLINE = System.getProperty("line.separator");
	
	private ArrayList<StatLoF> first, second;
	String outfile;
	StringBuilder result;
	
	public Comparison(String file1, String file2, String outfile) {
		try {
			first = readStatFile(file1);
			second = readStatFile(file2);
			this.outfile = outfile;
			
			String toolA = new File(file1).getName();
			String toolB = new File(file2).getName();
			
			result = new StringBuilder();
			
			//append header
			result.append("chr\tpos\tref_allele\talt_allele\t");
			result.append("effect_in_"+toolA+"\t");
			result.append("effect_in_"+toolB+"\t");
			result.append("consequence_in_"+toolA+"\t");
			result.append("consequence_in_"+toolB);
			result.append(NEWLINE);
			
		} catch (IOException e) {
			System.err.println(file1+ " could not be read!");
			e.printStackTrace();
		}
	}
	
	public void compare() {
		this.calcDiff();
		try {
			this.writeDiff();
		} catch (IOException e) {
			System.err.println(outfile+" could not be written.");
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
			
			result.add(new StatLoF(chr, pos, ref_allele, alt_allele, consequences, effects));
		}
		return result;
	}
	
	private void calcDiff() {
		int i = 0;
		int j = 0;
		int comp = -1;
		while (i < first.size() && j < second.size()) {
			StatLoF lof1 = first.get(i);
			StatLoF lof2 = second.get(j);
			comp = lof1.comparePos(lof2);

			if (comp == -1) {
				result.append(lof1.getChr()+"\t"+lof1.getPos()+"\t"+lof1.getRef_allele()+"\t"+lof1.getAlt_allele()+"\t");
				result.append(lof1.getEffect()+"\t"+"no_variant"+"\t");
				result.append(lof1.getConsequence()+"\t"+"no_variant");
				result.append(NEWLINE);
				i++;
			} else if (comp == 1) {
				result.append(lof2.getChr()+"\t"+lof2.getPos()+"\t"+lof2.getRef_allele()+"\t"+lof2.getAlt_allele()+"\t");
				result.append("no_variant"+"\t"+lof2.getEffect()+"\t");
				result.append("no_variant"+"\t"+lof2.getConsequence());
				result.append(NEWLINE);
				j++;
			} else {
				//create two lists different reference alleles at same position and compare each
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
						result.append(mylof1.getChr()+"\t"+mylof1.getPos()+"\t"+mylof1.getRef_allele()+"\t"+mylof1.getAlt_allele()+"\t");
						result.append(mylof1.getEffect()+"\t"+mylof2.getEffect()+"\t");
						result.append(mylof1.getConsequence()+"\t"+mylof2.getConsequence());
						result.append(NEWLINE);
					}
				}
			}
		}
		while(i < first.size()) {
			StatLoF lof = first.get(i);
			result.append(lof.getChr()+"\t"+lof.getPos()+"\t"+lof.getRef_allele()+"\t"+lof.getAlt_allele()+"\t");
			result.append(lof.getEffect()+"\t"+"no_variant"+"\t");
			result.append(lof.getConsequence()+"\t"+"no_variant");
			result.append(NEWLINE);
			i++;
		}
		while(j < second.size()) {
			StatLoF lof = first.get(j);
			result.append(lof.getChr()+"\t"+lof.getPos()+"\t"+lof.getRef_allele()+"\t"+lof.getAlt_allele()+"\t");
			result.append(lof.getEffect()+"\t"+"no_variant"+"\t");
			result.append(lof.getConsequence()+"\t"+"no_variant");
			result.append(NEWLINE);
			j++;
		}
	}
	
	private void writeDiff() throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		bw.write(result.toString());
		bw.close();
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String file1,file2,out;
		
		file1 = "/home/tim/Dropbox/Studium/Helmholtz/Master/MacArthur/phase1/chr19_macArthur_phase1_v3.vcf";
		file2 = "/home/tim/LOF_Project_local/files/chr19_20150414/LOFTEE/chr19_20101123_DEL_replaced.VEP_annotation.lof_statistic.txt";
		out = "/home/tim/LOF_Project_local/files/MacArthurVAT_filter_vs_LOFTEE_all_chr19_phase1_v3.csv";
		
//		file1 = "/home/ibis/tim.jeske/KORAdata/analysis_ready/20150319_vat_k/analysis_ready.diabetes.filtered.haplotypecaller.vat.lof_statistic.txt";
//		file2 = "/home/ibis/tim.jeske/KORAdata/analysis_ready/20150401_loftee_q/all_lof_stats/variant_effect_LoF.lof_statistic.txt";
		
		Comparison my = new Comparison(file1,file2,out);
		my.compare();

	}

}

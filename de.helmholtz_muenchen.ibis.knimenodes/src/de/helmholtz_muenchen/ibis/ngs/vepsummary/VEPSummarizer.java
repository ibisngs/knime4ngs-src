/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package de.helmholtz_muenchen.ibis.ngs.vepsummary;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import de.helmholtz_muenchen.ibis.utils.ngs.AnnotationParser;
import de.helmholtz_muenchen.ibis.utils.ngs.BioEntity;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFFile;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFVariant;

/**
* @author Tim Jeske
*/

public class VEPSummarizer {
	
	public static void getVarSum(VCFFile vcf, AnnotationParser ap, HashMap<String, Gene> genes, HashMap<String,Boolean> sampleid2case, String outfile) throws IOException {
		
		//generate set of cases and controls
		HashSet<String> cases = new HashSet<>();
		HashSet<String> ctrls = new HashSet<>();
		for(String s: sampleid2case.keySet()) {
			if(sampleid2case.get(s)) {
				cases.add(s);
			} else {
				ctrls.add(s);
			}
		}
		
		//initialize header
		String header = "chr\tpos\tid\tref_allele\talt_allele\tAF\tobs_het\tobs_hom\tgene_id\tgene_symbol\teffect\tlof_trans\tall_trans\tconsequence\tobs_hom_case\tobs_hom_ctrl\tExAC_NFE_AF";
		
		//initialize writer
		BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
		bw.write(header);
		bw.newLine();
		
		//iterate over variants and convert a VCFVariant to a line for summary
		VCFVariant var;
		String var_line,alt_line, anno, gene_ids, gene_symbols, effects, lof_trans, all_trans;
		String [] alt_alleles;
		HashMap<String, HashSet<String>> gene2transcripts;
		LinkedList<Annotation> annos = new LinkedList<>();
		double af;
		
		while(vcf.hasNext()) {
			var = vcf.next();
			var_line = var.getChrom()+"\t"+var.getPos()+"\t"+var.getId()+"\t"+var.getRef()+"\t";
			anno = var.getInfoField(ap.getAnnId());
			alt_alleles = var.getAlt().split(",");
			for(String id: ap.getAnnotatedAlleleIds(anno)) {
				int i = Integer.parseInt(id);
				af = var.getAF(i);
				if(af>0.0) {
					
					gene_ids = "";
					gene_symbols = "";
					effects = "";
					lof_trans = "";
					all_trans = "";
					
					gene2transcripts = ap.getGene2TranscriptIds(i, anno);
					for(String g: gene2transcripts.keySet()) {
						gene_ids+=","+g;
						if(genes.get(g)==null) continue;
						gene_symbols+=","+genes.get(g).getSymbol();
						lof_trans+=","+gene2transcripts.get(g).size();
						all_trans+=","+genes.get(g).getTranscripts().size();
						
						if(gene2transcripts.get(g).size() == genes.get(g).getTranscripts().size()) {
							effects+=",full";
						} else {
							effects+=",partial";
						}
					}
					
					gene_ids = gene_ids.replaceFirst(",", "");
					gene_symbols = gene_symbols.replaceFirst(",", "");
					effects = effects.replaceFirst(",", "");
					lof_trans = lof_trans.replaceFirst(",", "");
					all_trans = all_trans.replaceFirst(",", "");
					
					annos = ap.getAnnotations(i, anno);
					alt_line = alt_alleles[i-1]+"\t"+af+"\t"+var.getHetCount(i)+"\t"+var.getHomCount(i)+"\t"+gene_ids+"\t"+gene_symbols+"\t"+effects+"\t"+lof_trans+"\t"+all_trans+"\t"+annos.toString()+"\t"+var.getHomCount(i,cases)+"\t"+var.getHomCount(i,ctrls)+"\t"+ap.getExACAF(i, anno);
					bw.write(var_line + alt_line);
					bw.newLine();	
				}
			}
		}
		bw.close();
	}
	
//	public static void getGeneSum(VCFFile vcf, AnnotationParser ap, HashMap<String, Gene> genes, HashMap<String,Boolean> sampleid2case, String outfile) throws IOException {
//		
//		int nr_cases = 0;
//		int nr_controls = 0;
//		
//		for(boolean a: sampleid2case.values()) {
//			if(a) nr_cases++;
//			else nr_controls++;
//		}
//		
//		GeneSummary rs = new GeneSummary(vcf,ap, true);
//		HashMap<String, ContingencyTable> tables = rs.getTables(sampleid2case);
//		
//		//writing sample counts for each gene
//		
//		//initialize header
//		String header = "gene_id\tgene_symbol\taff_case\taff_ctrl\tun_case\tun_ctrl";
//
//		BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
//		bw.write(header);
//		bw.newLine();
//		
//		String line;
//		
//		
////		for(String gene: tables.keySet()) {
////			line = gene+"\t"+genes.get(gene).getSymbol()+"\t"+tables.get(gene).verticalToString();
////			bw.write(line);
////			bw.newLine();
////		}
//		
//		for(String gene: genes.keySet()) {
//			line = gene+"\t"+genes.get(gene).getSymbol();
//			if(tables.containsKey(gene)) {
//				line += "\t"+tables.get(gene).verticalToString();
//			} else {
//				line += "\t"+new ContingencyTable(0,nr_cases,0,nr_controls).verticalToString();
//			}
//			bw.write(line);
//			bw.newLine();
//		}
//		
//		bw.close();
//	}

	public static void getSampleSum(VCFFile vcf, AnnotationParser ap, HashMap<String, Gene> genes, String outfile) throws IOException {
		HashMap<String, HashMap<String, HashMap<String,String>>> sample2gene2transcript2aff = new HashMap<>();
		
		VCFVariant var;
		String anno;
		HashMap<String, HashSet<String>> gene2trans;
		HashMap<String, String> t2aff;
		HashMap<String, HashMap<String, String>> gene2t2aff;
		ArrayList<String> samples = vcf.getSampleIds();
		String aff = "undef";
		
		while(vcf.hasNext()) {
			var = vcf.next();
			anno = var.getInfoField(ap.getAnnId());
			
			for(String id: ap.getAnnotatedAlleleIds(anno)) {
				gene2trans = ap.getGene2TranscriptIds(Integer.parseInt(id), anno);
				
				for(String s: samples) {
					aff = var.getAff(s,Integer.parseInt(id));
					if(!aff.equals("undef") && !aff.equals("unaff")) {
//					if(var.isAffected(s,Integer.parseInt(id))) {
						
						if(sample2gene2transcript2aff.containsKey(s)) {
							gene2t2aff = sample2gene2transcript2aff.get(s);
						} else {
							gene2t2aff = new HashMap<>();
							sample2gene2transcript2aff.put(s, gene2t2aff);
						}
						
						for(String g: gene2trans.keySet()) {
							
							if(gene2t2aff.containsKey(g)) {
								t2aff = gene2t2aff.get(g);
							} else {
								t2aff = new HashMap<>();
								gene2t2aff.put(g, t2aff);
							}
							
							for(String t: gene2trans.get(g)) {
								if(t2aff.containsKey(t)) {
									t2aff.put(t, getAff(t2aff.get(t), aff));
								} else {
									t2aff.put(t,aff);
								}
							}
						}
					}
				}
			}
		}
		
		BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
		String geneanno,ko_genes;
		boolean is_ko;
		boolean oneInactive;
		int at_least_one, ko;
		bw.write("sample\tnr_aff\tat_least_one_transcript_complete\tnr_complete\tgenes_complete\taffected_genes");
		bw.newLine();
		
		for(String s: sample2gene2transcript2aff.keySet()) {
			gene2t2aff = sample2gene2transcript2aff.get(s);
			geneanno = "";
			ko_genes = "";
			at_least_one = 0;
			ko = 0;
			
			for(String g: gene2t2aff.keySet()) {
				t2aff = gene2t2aff.get(g);
				is_ko = true;
				oneInactive = false;
				
				for(String t: t2aff.keySet()) {
					geneanno += ";"+"("+g+","+t+","+t2aff.get(t)+")";
					if(t2aff.get(t).equals("hom") || t2aff.get(t).equals("comp")) {
						oneInactive=true;
					} else {
						is_ko = false;
					}
				}
					
				if((t2aff.keySet().size() == genes.get(g).getTranscripts().size()) && is_ko) { //all protein coding transcripts are affected
					ko++;
					ko_genes += ";"+genes.get(g).getSymbol();
				}
				
				if(oneInactive) {
					at_least_one++;
				}
			}
			
			geneanno = geneanno.replaceFirst(";", "");
			ko_genes = ko_genes.replaceFirst(";", "");
			bw.write(s+"\t"+gene2t2aff.keySet().size()+"\t"+at_least_one+"\t"+ko+"\t"+ko_genes+"\t"+geneanno);
			bw.newLine();
		}
		
		bw.close();
	}
	
	/*
	 * write output as matrix, samples as header, new row for each transcript
	 */
	public static void getMatrix(VCFFile vcf, AnnotationParser ap, HashMap<String, Gene> genes, HashMap<String,Boolean> sampleid2case,  String outfile) throws IOException {
		
		HashMap<Boolean,String> state = new HashMap<>();
		state.put(true, "case");
		state.put(false, "control");
		
		//create map transcript id -> gene id
		HashMap<String,String> transcript2gene = new HashMap<>();
		for(Gene g: genes.values()) {
			for(String t: g.getTranscripts()) {
				transcript2gene.put(t, g.getId());
			}
		}
		
		ArrayList<String> samples = vcf.getSampleIds();
		BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
		
		//write header
		StringBuilder header = new StringBuilder("transcript_id");
		for(String s:samples) {
			header.append("\t"+s+"_"+state.get(sampleid2case.get(s)));
		}
		bw.write(header.toString());
		bw.newLine();
		
		//create matrix
		HashMap<String, HashMap<String,String>> transcript2sample2aff = new HashMap<>();
		HashMap<String, String> sample2aff = null;
		
		VCFVariant var;
		String anno, aff, gene;
		String chr = "";
		StringBuilder tmp = new StringBuilder();
		
		while(vcf.hasNext()) {
			var = vcf.next();
			
			//write results
			if(!var.getChrom().equals(chr)) {
				chr = var.getChrom();
				
				for(String t: transcript2sample2aff.keySet()) {
					gene = transcript2gene.get(t);
					tmp = new StringBuilder(t+"_"+gene+"_"+genes.get(gene).getSymbol());
					sample2aff = transcript2sample2aff.get(t);
					for(String s: samples) {
						tmp.append("\t"+sample2aff.get(s));
					}
					bw.write(tmp.toString());
					bw.newLine();
				}
				transcript2sample2aff.clear();
			}
			
			anno = var.getInfoField(ap.getAnnId());
			
			HashMap<String, HashSet<Integer>> t2allele_ids = ap.getEntity2AlleleIds(anno, BioEntity.TRANSCRIPT_ID);
			
			for(String t: t2allele_ids.keySet()) {
				sample2aff = null;
				for(Integer id: t2allele_ids.get(t)) {
					
					if(var.getAF(id)>0.0) {
						
						//get sample2aff
						if(transcript2sample2aff.containsKey(t)) {
							sample2aff = transcript2sample2aff.get(t);
						} else {
							sample2aff = new HashMap<>();
							for(String s: samples) {
								sample2aff.put(s, "unaff"); //default state
							}
						}
						
						for(String s: samples) {
							aff = var.getAff(s,id);
							sample2aff.put(s, getAff(sample2aff.get(s), aff));
						}
					}
				}
				
				if(sample2aff != null) {
					transcript2sample2aff.put(t, sample2aff);
				}
			}
		}
		
		for(String t: transcript2sample2aff.keySet()) {
			gene = transcript2gene.get(t);
			tmp = new StringBuilder(t+"_"+gene+"_"+genes.get(gene).getSymbol());
			sample2aff = transcript2sample2aff.get(t);
			for(String s: samples) {
				tmp.append("\t"+sample2aff.get(s));
			}
			bw.write(tmp.toString());
			bw.newLine();
		}
		bw.close();
	}
	
	private static String getAff(String a, String b) {
		if(a.equals("het") && b.equals("het")) {
			return "pot_comp";
		} else if(a.equals(b)) {
			return a;
		} else if(a.equals("unaff")) {
			return b;
		} else if(b.equals("unaff")) {
			return a;
		}  else if(a.equals("hom") || b.equals("hom")) {
			return "hom";
		} else if(a.equals("comp") || b.equals("comp")) {
			return "comp";
		} else if (a.equals("pot_comp") || b.equals("pot_comp")) {
			return "pot_comp";
		} else if(a.equals("het") || b.equals("het")) {
			return "het";
		} else if((a.equals("mat") && b.equals("pat")) ||(a.equals("pat") && b.equals("mat"))) {
			return "comp";
		}
		return "undef";
	}
	
//	public static void getGeneSetSum(VCFFile vcf, AnnotationParser ap, HashMap<String, HashSet<String>> set2genes, HashMap<String,Boolean> sampleid2case, String outfile) throws IOException {
//		
//		//writing sample counts for each gene
//		
//		//initialize header
//		String header = "set\tsize\taff_case\taff_ctrl\tun_case\tun_ctrl\taff_genes";
//				
//		//initialize writer
//		BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
//		bw.write(header);
//		bw.newLine();
//		
//		String line;
//		
//		GeneSummary rs = new GeneSummary(vcf,ap,false);
//		rs.groupBy(set2genes);
//		HashMap<String, ContingencyTable> tables = rs.getSetTables(sampleid2case);
//		
//		for(String set: tables.keySet()) {
//			line = set+"\t"+set2genes.get(set).size()+"\t"+tables.get(set).verticalToString()+"\t"+set2genes.get(set).toString().replaceAll("\\[|\\]", "");
//			bw.write(line);
//			bw.newLine();
//		}
//		bw.close();
//	}

}

package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import de.helmholtz_muenchen.ibis.ngs.caseControlAnalyzer.ContingencyTable;
import de.helmholtz_muenchen.ibis.utils.ngs.AnnotationParser;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFFile;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFVariant;

public class LOFSummarizer {
	
	public static void getVarSum(VCFFile vcf, AnnotationParser ap, HashMap<String, Gene> genes, String outfile) throws IOException {
		
		//initialize header
		String header = "chr\tpos\tid\tref_allele\talt_allele\tAF\tobs_het\tobs_hom\tgene_id\tgene_symbol\teffect\tlof_trans\tall_trans\tconsequence";
		
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
					alt_line = alt_alleles[i-1]+"\t"+af+"\t"+var.getHetCount(i)+"\t"+var.getHomCount(i)+"\t"+gene_ids+"\t"+gene_symbols+"\t"+effects+"\t"+lof_trans+"\t"+all_trans+"\t"+annos.toString();
					bw.write(var_line + alt_line);
					bw.newLine();	
				}
			}
		}
		bw.close();
	}
	
	public static void getGeneSum(VCFFile vcf, AnnotationParser ap, HashMap<String, Gene> genes, HashMap<String,Boolean> sampleid2case, String outfile) throws IOException {
		
		RegionSummary rs = new RegionSummary(vcf,ap, true);
		HashMap<String, ContingencyTable> tables = rs.getTables(sampleid2case);
		
		//writing sample counts for each gene
		
		//initialize header
		String header = "gene_id\tgene_symbol\taff_case\taff_ctrl\tun_case\tun_ctrl";

		BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
		bw.write(header);
		bw.newLine();
		
		String line;
		
		for(String gene: tables.keySet()) {
			line = gene+"\t"+genes.get(gene).getSymbol()+"\t"+tables.get(gene).verticalToString();
			bw.write(line);
			bw.newLine();
		}
		bw.close();
	}

	public static void getSampleSum(VCFFile vcf, AnnotationParser ap, HashMap<String, Gene> genes, String outfile) throws IOException {
		HashMap<String, HashMap<String, HashMap<String,String>>> sample2gene2transcript2aff = new HashMap<>();
		
		VCFVariant var;
		String anno;
		HashMap<String, HashSet<String>> gene2trans;
		HashMap<String, String> t2aff;
		HashMap<String, HashMap<String, String>> gene2t2aff;
		ArrayList<String> samples = vcf.getSampleIds();
		String aff;
		
		while(vcf.hasNext()) {
			var = vcf.next();
			anno = var.getInfoField(ap.getAnnId());
			
			for(String id: ap.getAnnotatedAlleleIds(anno)) {
				gene2trans = ap.getGene2TranscriptIds(Integer.parseInt(id), anno);
				
				for(String s: samples) {
					if(!(aff = var.getAff(s,Integer.parseInt(id))).equals("undef")) {
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
		String geneanno, final_aff;
		int hom, ko;
		bw.write("sample\tnr_aff\tnr_hom\tnr_ko\taffected_genes");
		bw.newLine();
		
		for(String s: sample2gene2transcript2aff.keySet()) {
			gene2t2aff = sample2gene2transcript2aff.get(s);
			geneanno = "";
			hom = 0;
			ko = 0;
			
			for(String g: gene2t2aff.keySet()) {
				t2aff = gene2t2aff.get(g);
				final_aff = "";
				for(String t: t2aff.keySet()) {
					final_aff = getAff(final_aff, t2aff.get(t));
				}
				
				//counter for homozygously affected genes
				if(final_aff.equals("hom")) {
					hom++;
				}
				
				//check ko of gene
				if((final_aff.equals("hom") || final_aff.equals("comp")) && t2aff.keySet().size() == genes.get(g).getTranscripts().size()) {
					final_aff = "ko";
					ko++;
				}
				geneanno += ";"+"("+g+","+final_aff+")";
			}
			
			geneanno = geneanno.replaceFirst(";", "");
			bw.write(s+"\t"+gene2t2aff.keySet().size()+"\t"+hom+"\t"+ko+"\t"+geneanno);
			bw.newLine();
		}
		
		bw.close();
	}
	
	private static String getAff(String a, String b) {
		if(a.equals(b)) {
			return a;
		} else if(a.equals("")) {
			return b;
		} else if(b.equals("")) {
			return a;
		} else if(a.equals("hom") || b.equals("hom")) {
			return "hom";
		} else if(a.equals("comp") || b.equals("comp")) {
			return "comp";
		} else if(a.equals("het") || b.equals("het")) {
			return "het";
		} else if((a.equals("mat") && b.equals("pat")) ||(a.equals("pat") && b.equals("mat"))) {
			return "comp";
		}
		System.out.println("returned undef for: "+a+" "+b);
		return "undef";
	}
	
	public static void getEnrichmentSum(VCFFile vcf, AnnotationParser ap, HashMap<String, HashSet<String>> set2genes, HashMap<String,Boolean> sampleid2case, String outfile) throws IOException {
		
		//writing sample counts for each gene
		
		//initialize header
		String header = "set\taff_case\taff_ctrl\tun_case\tun_ctrl";
				
		//initialize writer
		BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
		bw.write(header);
		bw.newLine();
		
		String line;
		
		RegionSummary rs = new RegionSummary(vcf,ap,false);
		rs.groupBy(set2genes);
		HashMap<String, ContingencyTable> tables = rs.getSetTables(sampleid2case);
		
		for(String set: tables.keySet()) {
			line = set+"\t"+tables.get(set).verticalToString();
			bw.write(line);
			bw.newLine();
		}
		bw.close();
	}

}

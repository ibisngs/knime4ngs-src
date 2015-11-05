package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

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
		
		HashMap<String, HashSet<String>> unaffected_samples = new HashMap<>();
		HashMap<String, HashSet<String>> affected_samples = new HashMap<>();
		HashMap<String, HashSet<Integer>> gene2allele_ids;
		HashSet<Integer> alt_alleles;
		HashSet<String> aff_samp, unaff_samp;
		VCFVariant var;
		String anno;
		
		while(vcf.hasNext()) {
			var = vcf.next();
			anno = var.getInfoField(ap.getAnnId());
			gene2allele_ids = ap.getGene2AlleleIds(anno);
			
			for(String gene: gene2allele_ids.keySet()) {
				alt_alleles = gene2allele_ids.get(gene);
				aff_samp = var.getAffectedSamples(alt_alleles);
				unaff_samp = var.getUnaffectedSamples(alt_alleles);
				if(unaffected_samples.containsKey(gene)) {
					//unaffected samples contains all samples that carry only non-LOF alleles in a gene
					unaffected_samples.get(gene).retainAll(unaff_samp);
				} else {
					unaffected_samples.put(gene,unaff_samp);
				}
				
				if(affected_samples.containsKey(gene)) {
					affected_samples.get(gene).addAll(aff_samp);
				} else {
					affected_samples.put(gene, aff_samp);
				}
			}
		}
		
		//writing sample counts for each gene
		
		//initialize header
		String header = "gene_id\tgene_symbol\taff_case\taff_ctrl\tun_case\tun_ctrl";
				
		//initialize writer
		BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
		bw.write(header);
		bw.newLine();
		
		String line;
		int aff_case, aff_ctrl, un_case, un_ctrl;
		for(String gene: affected_samples.keySet()) {
			
			aff_case = 0;
			aff_ctrl = 0;
			un_case = 0;
			un_ctrl = 0;
			
			for(String s: affected_samples.get(gene)) {
				if(sampleid2case.get(s)) {
					aff_case++;
				} else {
					aff_ctrl++;
				}
			}
			
			for(String s: unaffected_samples.get(gene)) {
				if(sampleid2case.get(s)) {
					un_case++;
				} else {
					un_ctrl++;
				}
			}
			
			if(aff_case>0 || aff_ctrl>0) {
				line = gene+"\t"+genes.get(gene).getSymbol()+"\t"+aff_case+"\t"+aff_ctrl+"\t"+un_case+"\t"+un_ctrl;
				bw.write(line);
				bw.newLine();
			}
		}
		bw.close();
	}

	public static void getSampleSum(VCFFile vcf, AnnotationParser ap, String outfile) throws IOException {
		HashMap<String, HashMap<String, HashMap<String,String>>> sample2gene2transcript2aff = new HashMap<>();
		
		VCFVariant var;
		String anno;
		HashMap<String, HashSet<String>> gene2trans;
		HashMap<String, String> t2aff;
		HashMap<String, HashMap<String, String>> gene2t2aff;
		Set<String> samples;
		
		while(vcf.hasNext()) {
			var = vcf.next();
			samples = var.getSampleIds();
			anno = var.getInfoField(ap.getAnnId());
			
			for(String id: ap.getAnnotatedAlleleIds(anno)) {
				gene2trans = ap.getGene2TranscriptIds(Integer.parseInt(id), anno);
				
				for(String s: samples) {	
					if(var.isAffected(s,Integer.parseInt(id))) {
						
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
								t2aff.put(t,"aff");
							}
						}
					}
				}
			}
		}
		
		BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
		
		bw.write("sample\tnr_aff\taffected_genes");
		bw.newLine();
		
		for(String s: sample2gene2transcript2aff.keySet()) {
			bw.write(s+"\t"+sample2gene2transcript2aff.get(s).keySet().size()+"\t"+sample2gene2transcript2aff.get(s).keySet());
			bw.newLine();
		}
		
		bw.close();
	}
}

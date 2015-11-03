package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

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
	
	public static void getGeneSum(VCFFile vcf, AnnotationParser ap, String outfile) {
		
	}

	public static void getSampleSum(VCFFile vcf, AnnotationParser ap, String outfile) {
	
	}

}

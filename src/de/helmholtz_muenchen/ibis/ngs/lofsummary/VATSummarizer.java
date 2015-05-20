package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.ArrayList;
import java.util.HashMap;

public class VATSummarizer extends Summarizer{

	//map VAT terms to SO terms
	HashMap<String,String> VATtoSO;
	
	public VATSummarizer(String vcf_file, String cds_file) {
		super(vcf_file, cds_file);
		VATtoSO = new HashMap<>();
		VATtoSO.put("spliceOverlap", "splice_site_variant");
		VATtoSO.put("synonymous", "synonymous_variant");
		VATtoSO.put("endOverlap", "coding_sequence_variant");
		VATtoSO.put("svOverlap", "structural_variant");
		VATtoSO.put("removedStop", "stop_lost");
		VATtoSO.put("prematureStop", "stop_gained");
		VATtoSO.put("deletionNFS", "inframe_deletion");
		VATtoSO.put("insertionNFS", "inframe_insertion");
		VATtoSO.put("nonsynonymous", "missense_variant");
		VATtoSO.put("startOverlap", "initiator_codon_variant");
		VATtoSO.put("insertionFS", "frameshift_variant");
		VATtoSO.put("deletionFS", "frameshift_variant");
	}

	@Override
	void getLoFVariant(String chr, String pos, String ref, String id, String alt, String af, int GT_index, String[] infos, ArrayList<String> genotypes) {
		HashMap<String,LoFGene> genes = new HashMap<String,LoFGene>();
		
		int nrOfTitles = additional_titles.size();
		
		//init additional fields
		
		int observed_homo = 0;
		int observed_hetero = 0;
		String gene_id= "";
		String consequence="";
		String gene_symbol;
		
		String [] allele_annotations = null;
		String [] annotation_fields = null;
		
		for(String i: infos) {
			if(i.startsWith("VA")) {
				i = i.replaceFirst("VA=","");
				allele_annotations = i.split(",");
				break;
			}
		}
		
		if(allele_annotations == null) return; //nothing was annotated
		
		//AlleleNumber:GeneName:GeneId:Strand:Type:FractionOfTranscriptsAffected:{List of transcripts}
		
		for(String a: allele_annotations) {
			annotation_fields = a.split(":");
			if(annotation_fields[0].equals(GT_index+"")) {
				consequence = VATtoSO.get(annotation_fields[4]);
				if(consequence == null) {
					consequence = annotation_fields[4];
				}
				gene_id = annotation_fields[2];
				if(gene_id.contains(".")) {
					gene_id = gene_id.split("\\.")[0];
				}
				gene_symbol = annotation_fields[1];
				gene_id2gene_symbol.put(gene_id, gene_symbol);
				
				LoFGene g;
				if(genes.containsKey(gene_id)) {
					g = genes.get(gene_id);
				} else {
					g = new LoFGene(nrOfTitles);
				}
				g.addConsequence(consequence);
				
				for(String s : annotation_fields) {
					if(s.startsWith("ENST")) {
						if(s.contains(".")) {
							g.addTranscript(s.split("\\.")[0]);
						} else {
							g.addTranscript(s);
						}
					}
				}
				genes.put(gene_id, g);
			}
		}
		
		if(genes.size()==0) return; //no LoF genes have been found
		
		ArrayList<String> adaptedGenotypes = new ArrayList<>();
		char gt = String.valueOf(GT_index).charAt(0);
		for(int i = 0; i<genotypes.size(); i++) {
			String str = genotypes.get(i);
			if(str.charAt(0)==str.charAt(2) && str.charAt(0)==gt) {
				observed_homo++;
				adaptedGenotypes.add(i,"1"+str.charAt(1)+"1");
			} else if (str.charAt(0)==gt || str.charAt(2)==gt) {
				observed_hetero++;
				adaptedGenotypes.add(i,str.replace(gt, '1'));
			} else {
				adaptedGenotypes.add(i,str);
			}
		}
		lof_statistic.add(new LoFVariant(chr, pos, ref, alt, id, af, observed_homo, observed_hetero, genes, genotypes));
		
	}
}

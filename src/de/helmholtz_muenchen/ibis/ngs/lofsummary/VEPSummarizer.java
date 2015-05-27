package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.ArrayList;
import java.util.HashMap;

public class VEPSummarizer extends Summarizer{
	
	public VEPSummarizer(String vcf_file, String cds_file) {
		super(vcf_file, cds_file);
	}

	@Override
	void getLoFVariant(String chr, String pos, String ref, String id, String alt, int GT_index, String[] infos, ArrayList<String> genotypes) {
		HashMap<String,LoFGene> genes = new HashMap<String,LoFGene>();
		
		int nrOfTitles = additional_titles.size();
		
		//init additional fields
		
		int observed_homo = 0;
		int observed_hetero = 0;
		String [] transcript_annotations=null;
		String [] annotation_fields;
		String gene_id= "";
		String gene_symbol;
		String consequence, transcript_id;
		
		for(String i:infos) {
			if(i.startsWith("CSQ")) {
				i = i.replaceFirst("CSQ=","");
				transcript_annotations = i.split(",");
				break;
			}
		}
		
		if(transcript_annotations == null) return; //nothing was annotated
		
		//fill HashMap genes
		for(String t: transcript_annotations) {
			
			annotation_fields = t.split("\\|");
			
			//work around for annotations ending with ...||
			String [] tmp = new String [nrOfVEPfields];
			for(int i = 0; i< annotation_fields.length;i++) {
				tmp[i] = annotation_fields[i];
			}
			for(int j = annotation_fields.length; j < tmp.length; j++) {
				tmp[j] = "";
			}
			annotation_fields = tmp;
			
			//alt_allele - ref_allele = annotation_fields[0];
			int allele_index = vep_header.get("allele");
//			if(alt.contains(annotation_fields[allele_index]) || annotation_fields[allele_index].equals("-")) { //found a LoF annotation
			if(annotation_fields[allele_index].equals(GT_index+"")) {
				consequence = annotation_fields[vep_header.get("consequence")];
				gene_symbol = annotation_fields[vep_header.get("gene_symbol")];
				gene_id = annotation_fields[vep_header.get("gene_id")];
				if(gene_id == null || gene_id.equals("")) continue;
				transcript_id = annotation_fields[vep_header.get("transcript_id")];
				
				LoFGene g;
				if(genes.containsKey(gene_id)) {
					g = genes.get(gene_id);
				} else {
					g = new LoFGene(nrOfTitles);
				}
				g.addConsequence(consequence);
				g.addTranscript(transcript_id);
				
				for(String s: additional_titles) {
					g.addInfo(annotation_fields[vep_header.get(s)], additional_titles.indexOf(s));
				}
				
				genes.put(gene_id, g);
				gene_id2gene_symbol.put(gene_id, gene_symbol);
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
		lof_statistic.add(new LoFVariant(chr, pos, ref, alt, id, observed_homo, observed_hetero, genes, genotypes));
	}

}

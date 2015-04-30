package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class VEPSummarizer extends Summarizer{

	boolean loftee = false;
	boolean exac = false;
	boolean cadd = false;
	
	public VEPSummarizer(String vcf_file, String cds_file, HashSet<String> usedPlugins) {
		super(vcf_file, cds_file);
		
		if(usedPlugins.contains("LOFTEE")) {
			loftee = true;
			additional_titles.add("confidence");
			additional_titles.add("failed_filters");
			additional_titles.add("lof_flags");
			additional_titles.add("lof_info");
		}
		if(usedPlugins.contains("ExAC")) {
			exac = true;
			additional_titles.add("ExAC_AF");
		}
		if(usedPlugins.contains("CADD")) {
			cadd = true;
			additional_titles.add("CADD_PHRED");
		}
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
			if(alt.contains(annotation_fields[allele_index]) || annotation_fields[allele_index].equals("-")) { //found a LoF annotation
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
				
				if(loftee) {
					String confidence = annotation_fields[vep_header.get("confidence")];
					String filter = annotation_fields[vep_header.get("lof_filter")];
					String lof_info = annotation_fields[vep_header.get("lof_info")];
					String lof_flags = annotation_fields[vep_header.get("lof_flags")];
					g.addInfo(confidence, additional_titles.indexOf("confidence"));
					g.addInfo(filter, additional_titles.indexOf("failed_filters"));
					g.addInfo(lof_info, additional_titles.indexOf("lof_info"));
					g.addInfo(lof_flags, additional_titles.indexOf("lof_flags"));
				}
				
				if(exac) {
					String exac_score = annotation_fields[vep_header.get("ExAC_AF")];
					g.addInfo(exac_score, additional_titles.indexOf("ExAC_AF"));
				}
				
				if(cadd) {
					String cadd_score = annotation_fields[vep_header.get("CADD_PHRED")];
					g.addInfo(cadd_score, additional_titles.indexOf("CADD_PHRED"));
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

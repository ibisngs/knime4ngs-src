package de.helmholtz_muenchen.ibis.ngs.filtersummary.output;

import de.helmholtz_muenchen.ibis.ngs.filtersummary.genes.GeneInfo;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.interactions.ReactomePPI;

import java.util.Iterator;
import java.util.LinkedList;

import de.helmholtz_muenchen.ibis.ngs.filtersummary.pathways.WikiPathway;
import de.helmholtz_muenchen.ibis.utils.lofs.Writer_Output;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf.VCFLine;

public class VariantFile {
	
	/*
	 * table = 1 -> do not output original sample and info field
	 * table = 2 -> write original sample and info field to output
	 */
	
	public static void createVarFile(String vcfs, LinkedList<VCFLine> x, GeneInfo gi, ReactomePPI rp, WikiPathway wp, String path, boolean lof, int table){
		
		Writer_Output wo = new Writer_Output(path);
		
		String header = "Type\t"+"Variant Effect\t"+"Tool\t"+"Chromosome\t"+"Position\t"+"Length\t"+
						"Reference Sequence\t"+"Variant Sequence\t"+"DBSNP ID\t"+"GATK Quality Score\t"+"Genotype\t"+"Supporting Reads\t"+
						"Coverage\t"+"GeneIDs (Ensembl)\t"+"Gene Names\t"+"Affected Transcripts\t"+"All Transcripts";
		
		if(gi!=null){
			header+="\tOMIM IDs";
		}
		if(rp!=null){
			header+="\tNumber of interactions (Reactome)";
		}
		if(wp!=null){
			header+="\tPathways (Wikipathway)";
		}
		
		if(table==2){
			header+= "\tInfo"+"\tSample";
		}
		header+="\t"+vcfs;
		
		wo.writeFileln(header);
		
		Iterator<VCFLine> it = x.iterator();
		
		while (it.hasNext()){
			VCFLine l = it.next();
			
			String type=l.getVariantType();
			String loftype=l.getLofs(lof);
			String tool=l.getToolName();
			String chrom=l.getChrom();
			String position=Integer.toString(l.getPosition());
			String length = Integer.toString(l.getLength());
			String refseq = l.getRef();
			String altseq = l.getAlt();
			String dbsnp = l.getID();
			String quality= l.getQual();
			String genotype = l.getGenotype();
			String reads = Integer.toString(l.getSuppReads());
			String coverage = Integer.toString(l.getCov());
			String ensemblids= l.getGenes(lof, 0);
			String genenames=l.getGenes(lof, 1);
			String afftrans=l.getGenes(lof, 2);
			String trans=l.getGenes(lof, 3);

			String row=type+"\t"+loftype+"\t"+tool+"\t"+chrom+"\t"+position+"\t"+length+"\t"
						+refseq+"\t"+altseq+"\t"+dbsnp+"\t"+quality+"\t"+genotype+"\t"+reads+"\t"
						+coverage+"\t"+ensemblids+"\t"+genenames+"\t"+afftrans+"\t"+trans;
			
			if(gi!=null){
				String omimids = getOmim(gi, ensemblids);
				row+="\t"+omimids;
			}
			
			if(rp!=null){
				String interactions=getDegree(rp, ensemblids);
				row+="\t"+interactions;
			}
			
			if(wp!=null){
				String pathways =getPathway(wp, ensemblids);
				row+="\t"+pathways;
			}
			
			if(table==2){
				String info = l.getInfo();
				String sample = l.getSample();
				row+="\t"+info+"\t"+sample;
			}
			
			wo.writeFileln(row);
		}
		
		wo.closew();
	}
	
	private static String getOmim(GeneInfo gi, String genelist){
		
		String [] split = genelist.split(",");
		String [] res = new String [split.length];
		
		for (int i =0; i<split.length; i++){
			String omim = gi.getOmimID(split[i]);
			if(omim.equals("")){
				res[i]="-";
			}
			else if(omim.contains(",")){
				if(split.length==1){
					res[i]=omim;
				}
				else{
					res[i]="("+omim+")";
				}
			}
			else{
				res[i]=omim;
			}
		}
		
		String omimlist = res[0];
		for(int i=1; i<res.length; i++){
			omimlist+=","+res[i];
		}
		
		return omimlist;
	}
	
	private static String getDegree(ReactomePPI rp, String genelist){
		
		String [] split = genelist.split(",");
		String [] res = new String[split.length];
		
		for(int i = 0; i<split.length; i++){
			if(rp.getHashMapAll().containsKey(split[i])){
				res[i]=Integer.toString(rp.getHashMapAll().get(split[i]).size());
			}
			else{
				res[i]="-";
			}
		}
		
		String degreelist = res[0];
		for(int i=1; i<res.length; i++){
			degreelist+=","+res[i];
		}
		
		return degreelist;
	}
	
	private static String getPathway(WikiPathway wp, String genelist){
		
		String [] split = genelist.split(",");
		String [] res = new String[split.length];
		
		for(int i = 0; i<split.length; i++){
			String pw = wp.getPathway(split[i]);
			if(pw.equals("")){
				res[i]="-";
			}
			else if(pw.contains(",")){
				if(split.length==1){
					res[i]=pw;
				}
				else{
					res[i]="("+pw+")";
				}
			}
			else{
				res[i]=pw;
			}
		}
		
		String pathwaylist = res[0];
		for(int i=1; i<res.length; i++){
			pathwaylist+=","+res[i];
		}
		
		return pathwaylist;
	}

}

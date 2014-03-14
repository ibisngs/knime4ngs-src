package de.helmholtz_muenchen.ibis.ngs.filtersummary.output;

import de.helmholtz_muenchen.ibis.ngs.filtersummary.genes.GeneList;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.interactions.ReactomePPI;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;
import java.util.Map.Entry;

import de.helmholtz_muenchen.ibis.ngs.filtersummary.pathways.WikiPathway;
import de.helmholtz_muenchen.ibis.utils.lofs.Writer_Output;

public class FilesGenes {
	
	// write list of all ensembl ids of the genelist
	public static void writeGeneList (GeneList l, String f){
		Writer_Output w = new Writer_Output(f);
		
		Iterator<Entry<String, Integer>> i = l.geneCount.entrySet().iterator();
		while(i.hasNext()){
			Entry<String, Integer> e=i.next();
			w.writeFileln(e.getKey().split("\\.")[0]);
		}
		
		w.closew();
	}
	
	// how many lof variants are located in a gene (ensembl id)
	public static void writeLofsperGene(GeneList l, String f){
		
		Writer_Output w = new Writer_Output(f);
		w.writeFileln("GeneID\t#LoFs");
	
		Iterator<Entry<String, Integer>> i = l.geneCount.entrySet().iterator();
		while(i.hasNext()){
			Entry<String, Integer> e=i.next();
			w.writeFileln(e.getKey()+"\t"+e.getValue());
		}
		
		w.closew();
		
	}
	
	// table with gene counts
	public static void writeGeneCountTable(GeneList l, String f){
		
		Writer_Output w = new Writer_Output(f);
		w.writeFileln("#genes\tgatk\tpindel\tfull\tomim");
		w.writeFileln(l.genes+"\t"+l.gatk+"\t"+l.pindel+"\t"+l.full+"\t"+l.omim);
		w.closew();
		
	}
	
	
	/* write interaction degrees of lof genes and all genes in reactome
	 * int interaction: interaction type
	 * 0 -> all interactions
	 * 1 -> direct complex
	 * 2 -> direct reaction 
	 */
	public static void writeGeneDegree(GeneList l, ReactomePPI rp, String genedegfile, String alldegfile, int interactions){
		
		if(rp==null){
			return;
		}
		
		HashMap<String, LinkedList<String>> hm=null;
		if(interactions==0){
			hm = rp.getHashMapAll();
		}
		else if(interactions==1){
			hm = rp.getHashMapComplex();
		}
		else if(interactions==2){
			hm= rp.getHashMapReaction();
		}
		
		Writer_Output wo = new Writer_Output(alldegfile);
		Iterator<Entry<String, LinkedList<String>>> it =hm.entrySet().iterator();
		while(it.hasNext()){
			Entry<String, LinkedList<String>> e = it.next();
			wo.writeFileln(e.getKey()+"\t"+e.getValue().size());
		}
		wo.closew();
		
		Writer_Output w = new Writer_Output(genedegfile);
		Set<String> s = l.geneCount.keySet();
		Iterator<String> i = s.iterator();
		while(i.hasNext()){
			String geneid=i.next().split("\\.")[0];
			
			if(hm.containsKey(geneid)){
				w.writeFileln(geneid+"\t"+hm.get(geneid).size());
			}
		}
		w.closew();		
	}
	
	// write pathways associated with lof genes
	public static void writePathways(GeneList l, WikiPathway wp, String file){
		
		if(wp==null){
			return;
		}
		
		Writer_Output w = new Writer_Output(file);
		Set<String> s = l.geneCount.keySet();
		Iterator<String> i = s.iterator();
		while(i.hasNext()){
			String geneid=i.next().split("\\.")[0];	
			if(wp.getHashMap().containsKey(geneid)){
				Iterator<String> j =wp.getHashMap().get(geneid).iterator();
				while(j.hasNext()){
					String pathway = j.next();
					w.writeFileln(geneid+"\t"+pathway);
				}
			}
		}
		w.closew();
	}
}

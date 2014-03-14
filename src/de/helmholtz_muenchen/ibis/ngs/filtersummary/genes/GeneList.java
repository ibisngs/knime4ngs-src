package de.helmholtz_muenchen.ibis.ngs.filtersummary.genes;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;
import java.util.Set;

import de.helmholtz_muenchen.ibis.ngs.filtersummary.pathways.WikiPathway;
import de.helmholtz_muenchen.ibis.utils.lofs.Writer_Output;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf.VCFFile;
import de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf.VCFLine;

public class GeneList {
	
	
	// hash map counts how many variants a gene contains
	public HashMap<String, Integer> geneCount;
	// hash map for tools and geneid
	public HashMap<String, Integer> geneTool;
	// hash map for full/ partial lof
	public HashMap<String, Boolean> geneFull;
	
	
	public int genes;
	public int pindel;
	public int gatk;
	public int full;
	public int omim;
	
	public GeneList(LinkedList<VCFLine> variants){
		
		geneCount=new HashMap<String, Integer>(700*2);
		geneTool=new HashMap<String, Integer>(700*2);
		geneFull=new HashMap<String, Boolean>(700*2);
		
		Iterator<VCFLine> it = variants.iterator();
		
		VCFLine l;
		while(it.hasNext()){

			l=it.next();
			Iterator<Gene> j = l.getGenes().iterator();
			
			Gene g;
			while(j.hasNext()){
				
				g=j.next();
				String lof=g.getType();
				String geneID=g.getID();
				
				// only consider genes with lof mutation
				if(lof.equals("prematureStop") || lof.equals("spliceOverlap") || lof.equals("removedStop") || 
				lof.equals("deletionFS") || lof.equals("insertionFS") || lof.equals("multiExonHit")){
			
			
					if(geneCount.containsKey(geneID)){
						geneCount.put(geneID, geneCount.get(geneID)+1);
					}
					else{
						geneCount.put(geneID, 1);
					}
					
					
					int tool = l.getTool();					
					if(geneTool.containsKey(geneID)){
						int hashtool = geneTool.get(geneID);
						
						if(hashtool!=VCFFile.GATK+VCFFile.PINDEL && hashtool!=tool){							
							geneTool.put(geneID, VCFFile.GATK+VCFFile.PINDEL);
						}
					}
					else{
						geneTool.put(geneID, tool);
					}
					
					
					boolean full=false;
					if(g.getTranscripts()==g.getAffTranscripts()){
						full=true;
					}
					
					if(geneFull.containsKey(geneID)){
						if(full){
							geneFull.put(geneID, full);
						}
					}
					else{
						geneFull.put(geneID, full);
					}
			
				}
				
			}
		}
	}
	
	public void countGenes (GeneInfo gi){
		
		// all genes
		genes=geneCount.size();
		
		// tool genes
		Iterator<Entry<String, Integer>> i = geneTool.entrySet().iterator();
		while(i.hasNext()){
			Entry<String, Integer> e=i.next();
			
			if(e.getValue()==VCFFile.GATK){
				gatk++;
			}
			else if(e.getValue()==VCFFile.PINDEL){
				pindel++;
			}
			else if(e.getValue()==VCFFile.GATK+VCFFile.PINDEL){
				gatk++;
				pindel++;
			}
			else{
				System.out.println("Invalid tool value");
			}
		}
		
		// full lof genes
		Iterator<Entry<String, Boolean>> j = geneFull.entrySet().iterator();
		while(j.hasNext()){
			Entry<String, Boolean> e=j.next();
			if(e.getValue()){
				full++;
			}
		}
		
		if(gi!=null){
			Set<String> s = geneCount.keySet();
			Iterator<String> k = s.iterator();
			while(k.hasNext()){
				String id = k.next();
				
				//remove last digits of ensemble identifier
				String [] split = id.split("\\.");
				
				String q=split[0];
				if(gi.ensembltomim.containsKey(q)){
					omim++;
				}
			}
		}
		else{
			omim=0;
		}
	}
		
	public void writePathways(WikiPathway wp, String file){
		
		Writer_Output w = new Writer_Output(file);
		
		Set<String> s = geneCount.keySet();
		Iterator<String> i = s.iterator();
		while(i.hasNext()){
			String geneid=i.next().split("\\.")[0];
			//System.out.println(geneid);
			
			if(wp.getHashMap().containsKey(geneid)){
				Iterator<String> j =wp.getHashMap().get(geneid).iterator();
				while(j.hasNext()){
					String pathway = j.next();
					//System.out.println(geneid+"\t"+pathway);
					w.writeFileln(geneid+"\t"+pathway);
				}
			}
		}
		
		w.closew();
	}


}

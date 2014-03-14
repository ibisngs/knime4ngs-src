package de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf;

import java.util.Iterator;

public class VCFFileGatkIndel extends VCFFile {
	
	private int dbsnps=-1;
	
	public VCFFileGatkIndel (String filename){
		super(filename, GATK);
	}
	
	private void countdbsnps(){
		Iterator<VCFLine> it = entries.iterator();
		
		VCFLine l;
		while(it.hasNext()){
			l=it.next();
			
			if(!l.getID().equals(".")){
				dbsnps++;
			}
		}
				
	}
	
	public int getdbsnpcount(){
		
		if(dbsnps==-1){
			countdbsnps();
		}
		return dbsnps;
	}
	
	public int getType(VCFLine line){
		
		String alt=line.getAlt();
		String ref=line.getRef();
		
		String [] split = alt.split(",");
		
		if(split.length==1){
			if(alt.length()>ref.length()){
				return VCFLine.INS;
			}
			else if(alt.length()<ref.length()){
				return VCFLine.DEL;
			}
			else{
				System.out.println("reflen==altlen");
				return 0;
			}
		}
		else if(split.length==2){
			boolean ins=false;
			boolean del=false;
			
			for (int i=0; i<2; i++){
				if(split[i].length()<ref.length()){
					del=true;
				}
				else if(split[i].length()>ref.length()){
					ins=true;
				}
				else{
					System.out.println("reflen==altlen");
				}
			}
			
			if(ins&& del){
				return VCFLine.INDEL;
			}
			if(ins &&!del){
				return VCFLine.INS;
			}
			if(!ins && del){
				return VCFLine.DEL;
			}
		}

		return 0;
	}

}

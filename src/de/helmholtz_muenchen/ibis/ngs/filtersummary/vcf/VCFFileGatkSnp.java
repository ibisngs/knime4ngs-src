package de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf;

import java.util.Iterator;

public class VCFFileGatkSnp extends VCFFile {
	
	private int dbsnps=-1;
	
	public VCFFileGatkSnp(String filename){
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
		return VCFLine.SNP;
	}
	


}

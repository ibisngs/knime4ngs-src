package de.helmholtz_muenchen.ibis.ngs.filtersummary.filter;

import java.util.Iterator;
import java.util.LinkedList;

import de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf.VCFLine;

public class CovFilter {
	
	private int mincov;
	private double frac;
	
	public CovFilter(int mincov, double frac){
		this.mincov=mincov;
		this.frac=frac;
	}
	
	public CovFilter(){
		mincov=10;
		frac=0.2;
	}
	
	public String toString(){
		return "minimum coverage\t"+mincov+"\t fraction supporting reads\t"+frac;
	}
	
	public void setFilterOptions(int mincov, double frac){
		this.frac = frac;
		this.mincov=mincov;
	}

	public LinkedList<VCFLine> filter(LinkedList<VCFLine> v){
		
		LinkedList<VCFLine> res = new LinkedList<VCFLine>(); 
		
		Iterator<VCFLine> i = v.iterator();
		while(i.hasNext()){
			VCFLine l=i.next();
			
			if(l.getCov()>=mincov && ((double) l.getSuppReads())/l.getCov()>=frac){
				res.add(l);
			}
			
		}
		
		return res;
	}

}

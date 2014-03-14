package de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf;

public class VCFFilePindelDel extends VCFFile {
	
	public VCFFilePindelDel (String filename){
		super(filename, PINDEL);
	}
	
	public int getType(VCFLine line){
		return VCFLine.DEL;
	}

}

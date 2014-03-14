package de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf;

public class VCFFilePindelIns extends VCFFile {

	public VCFFilePindelIns (String filename){
		super(filename, PINDEL);
	}
	
	public int getType(VCFLine line){
		return VCFLine.INS;
	}
}

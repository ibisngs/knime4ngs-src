package de.helmholtz_muenchen.ibis.utils.ngs.frost;

//package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.IOException;
import java.util.ArrayList;

public class VCF_Comparison {

	protected ArrayList<TrioData> refData = new ArrayList<>();
	protected ArrayList<TrioData> toCompData = new ArrayList<>();
	
	private String refFile;
	private String toCompFile;
	
	public VCF_Comparison (String refFile, String toCompFile) {
		setRefFile(refFile);
		setToCompFile(toCompFile);
	}

	public String getRefFile() {
		return refFile;
	}

	public void setRefFile(String refFile) {
		this.refFile = refFile;
	}

	public String getToCompFile() {
		return toCompFile;
	}

	public void setToCompFile(String toCompFile) {
		this.toCompFile = toCompFile;
	}
	
	/**
	 * @throws IOException 
	 * 
	 * 
	 */
	public void getVcfInfo () throws IOException { // always two files, ref and toComp
		if (VCF_DataCollection.checkContig(this.refFile, this.toCompFile)) {
			this.refData = VCF_DataCollection.extractVcf(this.refFile);
			this.toCompData = VCF_DataCollection.extractVcf(this.toCompFile);
		}
			
	}
}

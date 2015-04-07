package de.helmholtz_muenchen.ibis.ngs.lofstatistics;

import java.util.ArrayList;

public class SampleStat {
	
	int fullLOFs, partLOFs;
	ArrayList<String> complete_LOF_genes;
	ArrayList<String> part_LOF_genes;
	
	public int getFullLOFs() {
		return fullLOFs;
	}

	public int getPartLOFs() {
		return partLOFs;
	}

	public ArrayList<String> getComplete_LOF_genes() {
		return complete_LOF_genes;
	}

	public ArrayList<String> getPart_LOF_genes() {
		return part_LOF_genes;
	}

	public SampleStat(int fullLOFs, int partLOFs,ArrayList<String> complete_LOF_genes, ArrayList<String> part_LOF_genes) {
		this.fullLOFs = fullLOFs;
		this.partLOFs = partLOFs;
		this.complete_LOF_genes = complete_LOF_genes;
		this.part_LOF_genes = part_LOF_genes;
	}
	
	

}

package de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf;

import java.util.HashMap;

public class Chromosome implements Comparable<Chromosome>{
	
	private static String [] chrom_names={"chrM", "chr1","chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11","chr12",
			"chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"};
	
	private static Integer [] chrom_id={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
	
	private static HashMap<String, Integer> chroms = new HashMap<String, Integer>();
	
	static{
		for (int i= 0; i<chrom_names.length; i++){
			chroms.put(chrom_names[i], chrom_id[i]);
		}
	}
	
	private String chrom;
	
	public Chromosome(String x){
		this.chrom=x;
	}

	@Override
	public int compareTo(Chromosome o) {
		
		return chroms.get(chrom).compareTo(chroms.get(o.chrom));
	}
	
}

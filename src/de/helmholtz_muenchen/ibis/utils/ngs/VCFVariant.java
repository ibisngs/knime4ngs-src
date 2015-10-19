package de.helmholtz_muenchen.ibis.utils.ngs;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

public class VCFVariant {
	
	private String chrom, pos, id, ref, alt, qual, filter, info, format;
	private HashMap<String, String> sample_genotype;
	
	public VCFVariant() {
		sample_genotype = new HashMap<>();
	}
	
	public String getChrom() {
		return chrom;
	}

	public void setChrom(String chrom) {
		this.chrom = chrom;
	}

	public String getPos() {
		return pos;
	}

	public void setPos(String pos) {
		this.pos = pos;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getRef() {
		return ref;
	}

	public void setRef(String ref) {
		this.ref = ref;
	}

	public String getAlt() {
		return alt;
	}

	public void setAlt(String alt) {
		this.alt = alt;
	}

	public String getQual() {
		return qual;
	}

	public void setQual(String qual) {
		this.qual = qual;
	}

	public String getFilter() {
		return filter;
	}

	public void setFilter(String filter) {
		this.filter = filter;
	}

	public String getInfo() {
		return info;
	}

	public void setInfo(String info) {
		this.info = info;
	}

	public String getFormat() {
		return format;
	}

	public void setFormat(String format) {
		this.format = format;
	}
	
	public void addGenotype(String sample, String genotype) {
		this.sample_genotype.put(sample, genotype);
	}
	
	public HashSet<String> getUnaffectedSamples(Collection<Integer> allele_nums) {
		String gt;
		HashSet<String> result = new HashSet<>();
		for(String s: sample_genotype.keySet()) {
			gt = sample_genotype.get(s).split(":")[0];
			if (gt.charAt(0) != '.'
					&& !allele_nums.contains(Character.getNumericValue(gt.charAt(0)))
					&& gt.charAt(2) != '.'
					&& !allele_nums.contains(Character.getNumericValue(gt.charAt(2)))) {
				result.add(s);
			}
		}
		return result;
	}
	
	public HashSet<String> getAffectedSamples(Collection<Integer> allele_nums) {
		String gt;
		HashSet<String> result = new HashSet<>();
		for(String s: sample_genotype.keySet()) {
			gt = sample_genotype.get(s).split(":")[0];
			if(gt.length()!=3) {
				continue;
			}
			for(int i: allele_nums) {
				if((Character.getNumericValue(gt.charAt(0)) == i) || (Character.getNumericValue(gt.charAt(2)) == i)) {
					result.add(s);
				}
			}
		}
		return result;
	}
	
	public boolean isMultiallelic() {
		return alt.contains(",");
	}
	
	public String getInfoField(String field_id) {
		String [] fields = info.split(";");
		for(String f: fields) {
			if(f.startsWith(field_id)) {
				return f.split("=")[1];
			}
		}
		return null;
	}
}
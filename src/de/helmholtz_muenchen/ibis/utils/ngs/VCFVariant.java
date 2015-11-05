package de.helmholtz_muenchen.ibis.utils.ngs;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class VCFVariant {
	
	private String chrom, pos, id, ref, alt, qual, filter, info, format;
	private HashMap<String, String> sample_genotype;
	
	public VCFVariant() {
		sample_genotype = new HashMap<>();
	}
	
	public String getChrom() {
		return chrom.replaceAll("chr", "");
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
			if(gt.length()!=3) {
				continue;
			}
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
	
	public double getAF(int allele_id) {
		String gt;
		double an = 0.0;
		double ac = 0.0;
		for(String s: sample_genotype.values()) {
			gt = s.split(":")[0];
			if(gt.length()!=3) {
				continue;
			}
			if(gt.charAt(0)!='.') {
				an++;
				if(Character.getNumericValue(gt.charAt(0))==allele_id) {
					ac++;
				}
			}
			if(gt.charAt(2)!='.') {
				an++;
				if(Character.getNumericValue(gt.charAt(2))==allele_id) {
					ac++;
				}
			}
		}
		return ac/an;
	}
	
	public int getHomCount(int allele_id) {
		int res = 0;
		String gt;
		for(String s: sample_genotype.values()) {
			gt = s.split(":")[0];
			if(gt.length()!=3) {
				continue;
			}
			if(Character.getNumericValue(gt.charAt(0)) == allele_id && Character.getNumericValue(gt.charAt(2))== allele_id) {
				res++;
			}
		}
		return res;
	}
	
	public int getHetCount(int allele_id) {
		int res = 0;
		String gt;
		for(String s: sample_genotype.values()) {
			gt = s.split(":")[0];
			if(gt.length()!=3) {
				continue;
			}
			if(Character.getNumericValue(gt.charAt(0)) == allele_id ^ Character.getNumericValue(gt.charAt(2))== allele_id) {
				res++;
			}
		}
		return res;
	}
	
	public Set<String> getSampleIds()  {
		return this.sample_genotype.keySet();
	}
	
	public boolean isAffected(String sample, int id) {
		return sample_genotype.get(sample).split(":")[0].contains(id+"");
	}
}

package de.helmholtz_muenchen.ibis.ngs.annotationComparator;

import java.util.HashSet;

public class StatLoF {
	
	String chr, pos, ref_allele, alt_allele;
	HashSet<String> consequences, effects;
	
	public StatLoF(String chr, String pos, String ref_allele,
			String alt_allele, String[] consequences, String[] effects) {
		this.chr = chr;
		this.pos = pos;
		this.ref_allele = ref_allele;
		this.alt_allele = alt_allele;
		this.consequences = new HashSet<>();
		for(String s: consequences) {
			if(s.equals("splice_donor_variant") || s.equals("splice_acceptor_variant")) {
				s = "splice_site_variant";
			}
			this.consequences.add(s);
		}
		this.effects = new HashSet<>();
		for(String s: effects) {
			this.effects.add(s);
		}
	}

	public String getChr() {
		return chr;
	}

	public String getPos() {
		return pos;
	}

	public String getRef_allele() {
		return ref_allele;
	}

	public String getAlt_allele() {
		return alt_allele;
	}

	public HashSet<String> getConsequences() {
		return consequences;
	}

	public HashSet<String> getEffects() {
		return effects;
	}

	public int comparePos(StatLoF lof) {
		int a,b;
		if(!this.chr.equals(lof.getChr())) {
			try {
				a = Integer.parseInt(this.chr);
			} catch (NumberFormatException e) {
				return 1;
			}
			try {
				b = Integer.parseInt(lof.getChr());
			} catch (NumberFormatException e) {
				return -1;
			}
			if(a>b) return 1;
			if(a<b) return -1;
		}
		
		//chromosomes are equal
		a = Integer.parseInt(this.pos);
		b = Integer.parseInt(lof.pos);
		if(a>b) {
			return 1;
		} else if(a<b) {
			return -1;
		} else {
			return 0;
		}
	}
	
	public String toString() {
		return chr+"\t"+pos+"\t"+ref_allele+"\t"+alt_allele+"\t"+consequences+"\t"+effects;
	}
	
}

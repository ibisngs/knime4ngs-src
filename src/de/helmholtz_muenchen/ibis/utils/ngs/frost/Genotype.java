package de.helmholtz_muenchen.ibis.utils.ngs.frost;

/**
 * @author tanzeem.haque
 * Object for Mutation
 */

public class Genotype {

	private int allele1;
	private int allele2;
	
	
	public Genotype(int a1, int a2) {
		// TODO Auto-generated constructor stub
		setAllele1(a1);
		setAllele2(a2);
	}

	/**
	 * @return the allele1
	 */
	protected int getAllele1() {
		return allele1;
	}

	/**
	 * @param allele1 the allele1 to set
	 */
	protected void setAllele1(int allele1) {
		this.allele1 = allele1;
	}

	/**
	 * @return the allele2
	 */
	protected int getAllele2() {
		return allele2;
	}

	/**
	 * @param allele2 the allele2 to set
	 */
	protected void setAllele2(int allele2) {
		this.allele2 = allele2;
	}

}

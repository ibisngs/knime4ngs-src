package de.helmholtz_muenchen.ibis.utils.ngs.frost;
//package de.helmholtz_muenchen.ibis.utils.ngs.frost;


/**
 * @author tanzeem.haque
 * Object for Mutation
 *
 */
public class TrioData {
	
	private int position;
	private Allelic_Change alleles;
	private TrioGenotype trio;
	
	public TrioData(int position, Allelic_Change alleles, TrioGenotype trio) {
		// TODO Auto-generated constructor stub
		setPosition(position);
		setAlleles(alleles);
		setTrio(trio);
	}

	/**
	 * @return the position
	 */
	protected int getPosition() {
		return position;
	}

	/**
	 * @param position the position to set
	 */
	protected void setPosition(int position) {
		this.position = position;
	}

	/**
	 * @return the alleles
	 */
	protected Allelic_Change getAlleles() {
		return alleles;
	}

	/**
	 * @param alleles the alleles to set
	 */
	protected void setAlleles(Allelic_Change alleles) {
		this.alleles = alleles;
	}

	/**
	 * @return the parents
	 */
	protected TrioGenotype getTrio() {
		return trio;
	}

	/**
	 * @param parents the parents to set
	 */
	protected void setTrio(TrioGenotype trio) {
		this.trio = trio;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		// TODO Auto-generated method stub
		String data = "";
		data += getPosition() + "\t"
					+ getAlleles().getReference() + "\t"
					+ getAlleles().getAltered() + "\t"
					+ getTrio().getMother().getAllele1()
					+ "/"
					+ getTrio().getMother().getAllele2()
					+ "\t"
					+ getTrio().getFather().getAllele1()
					+ "/"
					+ getTrio().getFather().getAllele2()
					+ "\t"
					+ getTrio().getChild().getAllele1()
					+ "/"
					+ getTrio().getChild().getAllele2()
					+ "\n";
		
		return data;
	}

	
}

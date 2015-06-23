package de.helmholtz_muenchen.ibis.utils.ngs.frost;


/**
 * @author tanzeem.haque
 * Object for Mutation
 *
 */
public class ParentalData {
	
	private int position;
	private Allelic_Change alleles;
	private ParentalGenotype parents;
	private String parent;
	
	public ParentalData(int position, Allelic_Change alleles, ParentalGenotype parents, String parent) {
		// TODO Auto-generated constructor stub
		setPosition(position);
		setAlleles(alleles);
		setParents(parents);
		setParent(parent);
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
	protected ParentalGenotype getParents() {
		return parents;
	}

	/**
	 * @param parents the parents to set
	 */
	protected void setParents(ParentalGenotype parents) {
		this.parents = parents;
	}
	/**
	 * @return the parent
	 */
	public String getParent() {
		return parent;
	}

	/**
	 * @param parent the parent to set
	 */
	public void setParent(String parent) {
		this.parent = parent;
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
					+ getParents().getFather().getAllele1()
					+ "/"
					+ getParents().getFather().getAllele2()
					+ "\t"
					+ getParents().getMother().getAllele1()
					+ "/"
					+ getParents().getMother().getAllele2()
					+ "\t"
					+ getParent()
					+ "\n";
		
		return data;
	}

	
}

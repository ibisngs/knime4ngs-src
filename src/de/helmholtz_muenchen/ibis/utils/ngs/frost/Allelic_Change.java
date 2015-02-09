package de.helmholtz_muenchen.ibis.utils.ngs.frost;

/**
 * @author tanzeem.haque
 * Object for Mutation
 *
 */

public class Allelic_Change {

	private String reference;
	private String altered;
	
	public Allelic_Change(String reference, String altered) {
		// TODO Auto-generated constructor stub
		setReference(reference);
		setAltered(altered);
		
	}

	/**
	 * @return the reference
	 */
	protected String getReference() {
		return reference;
	}

	/**
	 * @param reference the reference to set
	 */
	protected void setReference(String reference) {
		this.reference = reference;
	}

	/**
	 * @return the altered
	 */
	protected String getAltered() {
		return altered;
	}

	/**
	 * @param altered the altered to set
	 */
	protected void setAltered(String altered) {
		this.altered = altered;
	}

}

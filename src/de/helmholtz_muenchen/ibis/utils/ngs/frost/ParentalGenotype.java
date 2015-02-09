package de.helmholtz_muenchen.ibis.utils.ngs.frost;


/**
 * @author tanzeem.haque
 *
 */
public class ParentalGenotype {

	private Genotype father;
	private Genotype mother;
	
	public ParentalGenotype(Genotype father, Genotype mother) {
		// TODO Auto-generated constructor stub
		setFather(father);
		setMother(mother);
	}

	/**
	 * @return the father
	 */
	public Genotype getFather() {
		return father;
	}

	/**
	 * @param father the father to set
	 */
	public void setFather(Genotype father) {
		this.father = father;
	}

	/**
	 * @return the mother
	 */
	public Genotype getMother() {
		return mother;
	}

	/**
	 * @param mother the mother to set
	 */
	public void setMother(Genotype mother) {
		this.mother = mother;
	}

}

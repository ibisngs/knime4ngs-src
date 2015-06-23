package de.helmholtz_muenchen.ibis.utils.ngs.frost;


/**
 * @author tanzeem.haque
 *
 */
public class ParentalGenotype {

	private Genotype mother;
	private Genotype father;

	
	public ParentalGenotype(Genotype mother, Genotype father) {
		// TODO Auto-generated constructor stub
		setMother(mother);
		setFather(father);
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

}

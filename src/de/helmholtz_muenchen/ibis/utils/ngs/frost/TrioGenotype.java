package de.helmholtz_muenchen.ibis.utils.ngs.frost;

//package de.helmholtz_muenchen.ibis.utils.ngs.frost;

public class TrioGenotype extends ParentalGenotype {
	
	private Genotype child;
	
	public TrioGenotype (Genotype father, Genotype mother, Genotype child) {
		super(father, mother);
		setChild(child);
	}

	public Genotype getChild() {
		return child;
	}

	public void setChild(Genotype child) {
		this.child = child;
	}

}
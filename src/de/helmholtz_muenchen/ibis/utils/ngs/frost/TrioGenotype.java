package de.helmholtz_muenchen.ibis.utils.ngs.frost;

//package de.helmholtz_muenchen.ibis.utils.ngs.frost;

public class TrioGenotype extends ParentalGenotype {
	
	private Genotype child;
	
	public TrioGenotype (Genotype mother, Genotype father, Genotype child) {
		super(mother, father);
		setChild(child);
	}

	public Genotype getChild() {
		return child;
	}

	public void setChild(Genotype child) {
		this.child = child;
	}

}
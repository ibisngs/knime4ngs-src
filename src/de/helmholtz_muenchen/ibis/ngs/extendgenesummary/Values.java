package de.helmholtz_muenchen.ibis.ngs.extendgenesummary;

public class Values {
	
	int lof_count;
	double p_lof;
	int max_degree;
	double max_cc;
	
	public Values(int lof_count, double p_lof, int max_degree, double max_cc) {
		super();
		this.lof_count = lof_count;
		this.p_lof = p_lof;
		this.max_degree = max_degree;
		this.max_cc = max_cc;
	}
	
	public int getLof_count() {
		return lof_count;
	}
	public void setLof_count(int lof_count) {
		this.lof_count = lof_count;
	}
	public double getP_lof() {
		return p_lof;
	}
	public void setP_lof(double p_lof) {
		this.p_lof = p_lof;
	}
	public int getMax_degree() {
		return max_degree;
	}
	public void setMax_degree(int max_degree) {
		if(max_degree > this.max_degree) {
			this.max_degree = max_degree;
		}
	}
	
	public double getMax_cc() {
		return max_cc;
	}
	public void setMax_cc(double max_cc) {
		if(max_cc > this.max_cc) {
			this.max_cc = max_cc;
		}
	}
	

}

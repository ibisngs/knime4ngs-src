package de.helmholtz_muenchen.ibis.ngs.caseControlAnalyzer;

public class ContingencyTable {
	
	//group		cond	ncond
	//case		a  		b
	//control	c  		d
	
	private int a,b,c,d;
	
	public ContingencyTable(int a, int b, int c, int d) {
		this.a = a;
		this.b = b;
		this.c = c;
		this.d = d;
	}

	public int getA() {
		return a;
	}

	public int getB() {
		return b;
	}

	public int getC() {
		return c;
	}

	public int getD() {
		return d;
	}

	public int getAll() {
		return a+b+c+d;
	}

}

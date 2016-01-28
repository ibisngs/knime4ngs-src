package de.helmholtz_muenchen.ibis.ngs.geneBasedAnalysis;

import java.util.ArrayList;
import java.util.Arrays;

public class CaseControlArray {

	short [] cases, controls;
	
	public CaseControlArray(short [] cases, short [] controls) {
		this.cases = cases;
		this.controls = controls;
	}
	
	public CaseControlArray(ArrayList<Short> cases, ArrayList<Short> controls) {
		this.cases = new short[cases.size()];
		for(int i = 0; i < this.cases.length; i++) {
			this.cases[i] = cases.get(i);
		}
		this.controls  = new short[controls.size()];
		for(int i = 0; i < this.controls.length; i++) {
			this.controls[i] = controls.get(i);
		}
	}
	
	public String toString() {
		String rep;
		
		rep = "cases: "+Arrays.toString(cases)+ System.getProperty("line.separator");
		rep += "controls: "+Arrays.toString(controls);
		
		return rep;
	}
	
	public short [] getCaseArray() {
		return this.cases;
	}
	
	public short [] getControlArray() {
		return this.controls;
	}
}

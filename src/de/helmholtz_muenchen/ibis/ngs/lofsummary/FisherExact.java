package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import org.apache.commons.math3.util.CombinatoricsUtils;

public class FisherExact {

    /**
     * calculates the P-value for this specific state
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return the P-value
     */
    public static double getP(int a, int b, int c, int d) {
        int n = a + b + c + d;
        
        double p = CombinatoricsUtils.binomialCoefficientDouble(a+b,b)/CombinatoricsUtils.binomialCoefficientDouble(n, a+c);
        p = p * CombinatoricsUtils.binomialCoefficientDouble(c+d, c);
        return p;
        
    }
    
    public static double oneTailedGreater(int a, int b, int c, int d) {
    	double p = 0.0;
    	int z = Math.min(b, c);
    	for(int i = 0; i <= z; i++) {
    		p += getP(a+i, b-i, c-i, d+i);
    	}
    	return p;
    }
    public static void main (String [] args) {
    	System.out.println(getP(1,9,11,3));
    	System.out.println(oneTailedGreater(1,9,11,3));
    	System.out.println(getP(104,105,0,0));
    	System.out.println(oneTailedGreater(104,105,0,0));
    }
}

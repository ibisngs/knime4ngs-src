package de.helmholtz_muenchen.ibis.ngs.caseControlAnalyzer;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.util.CombinatoricsUtils;

public class Statistics {

    /**
     * calculates the P-value for this specific state
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return the P-value
     */
    public static double getFisherP(int a, int b, int c, int d) {
        int n = a + b + c + d;
        
        double p = CombinatoricsUtils.binomialCoefficientDouble(a+b,b)/CombinatoricsUtils.binomialCoefficientDouble(n, a+c);
        p = p * CombinatoricsUtils.binomialCoefficientDouble(c+d, c);
        return p;
        
    }
    
    public static double getFisherOneTailedGreater(ContingencyTable ct) {
    	return getFisherOneTailedGreater(ct.getA(), ct.getB(), ct.getC(), ct.getD());
    }
    
    public static double getFisherOneTailedGreater(int a, int b, int c, int d) {
    	double p = 0.0;
    	int z = Math.min(b, c);
    	for(int i = 0; i <= z; i++) {
    		p += getFisherP(a+i, b-i, c-i, d+i);
    	}
    	return p;
    }
    
    public static double getBinomialBackground(ContingencyTable ct, double bg_freq) {
    	
    	int case_aff = ct.getA();
    	int n_cases = ct.getA()+ct.getB();
    	int control_aff = ct.getC();
    	int n_controls = ct.getC() + ct.getD();
    	
    	NormalDistribution nd = new NormalDistribution();
				
		double p_val_case_vs_bg = 1 - new BinomialDistribution(n_cases+1, bg_freq).cumulativeProbability(case_aff);
		double z_score_case_vs_bg = nd.inverseCumulativeProbability(p_val_case_vs_bg);
		
		double p_val_control_vs_bg = 1 - new BinomialDistribution(n_controls+1, bg_freq).cumulativeProbability(control_aff);
		double z_score_control_vs_bg = nd.inverseCumulativeProbability(p_val_control_vs_bg);

		double z_score_diff = z_score_case_vs_bg - z_score_control_vs_bg;

		return nd.cumulativeProbability(z_score_diff);
    }
    
    public static void main (String [] args) {
    	System.out.println(getFisherP(1,9,11,3));
    	System.out.println(getFisherOneTailedGreater(1,9,11,3));
    	System.out.println(getFisherP(104,105,0,0));
    	System.out.println(getFisherOneTailedGreater(104,105,0,0));
    }
}
package de.helmholtz_muenchen.ibis.ngs.caseControlAnalyzer;

import rcaller.RCaller;
import rcaller.RCode;

public class Statistics {
	
	private RCaller r = null;
	private RCode code = null;
	
	public Statistics() {
		r = new RCaller();
		r.setRscriptExecutable("/home/software/bin/Rscript");
		r.setRExecutable("/home/software/bin/R");
		code = new RCode();
		r.setRCode(code);
	}
	
	public void quit() {
		r.StopRCallerOnline();
		r.stopStreamConsumers();
	}

    /**
     * calculates the P-value for this specific state
     *
     * @param a     a, b, c, d are the four cells in a 2x2 matrix
     * @param b
     * @param c
     * @param d
     * @return the P-value
     */
//    public static double getFisherP(int a, int b, int c, int d) {
//        int n = a + b + c + d;
//        
//        double p = CombinatoricsUtils.binomialCoefficientDouble(a+b,b)/CombinatoricsUtils.binomialCoefficientDouble(n, a+c);
//        p = p * CombinatoricsUtils.binomialCoefficientDouble(c+d, c);
//        return p;
//        
//    }
    
    public double getFisherOneTailedGreater(ContingencyTable ct) {
    	
    	//add pseudo counts
    	int x = ct.getA()+1;
    	int m = x + ct.getB();
    	int n = ct.getC() + ct.getD();
    	int k = x + ct.getC();
    	
    	code.clear();
    	code.addRCode("res<-phyper("+(x-1)+","+m+","+n+","+k+",lower.tail=F)");
    	r.runAndReturnResultOnline("res");
    	double res = r.getParser().getAsDoubleArray("res")[0];
		r.deleteTempFiles();
    	return res;
    }
    
    public double getHypergeometricBackground(ContingencyTable ct, int pop_size, int pop_cases) {
    	int x_case = ct.getA()+1;
    	int x_ctrl = ct.getC()+1;
    	
    	int m_case = x_case + ct.getB();
    	int m_ctrl = x_ctrl + ct.getD();
    	
    	double pop_frac = (double)pop_cases/(double)pop_size;
    	double case_frac = (double)x_case/(double)m_case;
    	double ctrl_frac = (double)x_ctrl/(double)m_ctrl;
    	
    	code.clear();
    	
    	if(pop_frac > case_frac) {
        	code.addRCode("z_case_bg<--qnorm(phyper("+(x_case-1)+","+m_case+","+pop_size+","+(pop_cases+x_case)+",lower.tail=T))");
    	} else {
        	code.addRCode("z_case_bg<-qnorm(phyper("+(x_case-1)+","+m_case+","+pop_size+","+(pop_cases+x_case)+",lower.tail=F))");
    	}
    	
    	if(pop_frac > ctrl_frac) {
        	code.addRCode("z_ctrl_bg<--qnorm(phyper("+(x_ctrl-1)+","+m_ctrl+","+pop_size+","+(pop_cases+x_ctrl)+",lower.tail=T))");
    	} else {
        	code.addRCode("z_ctrl_bg<-qnorm(phyper("+(x_ctrl-1)+","+m_ctrl+","+pop_size+","+(pop_cases+x_ctrl)+",lower.tail=F))");
    	}
    	code.addRCode("res<-pnorm(z_case_bg-z_ctrl_bg)");

		r.runAndReturnResultOnline("res");
		double res = r.getParser().getAsDoubleArray("res")[0];
		r.deleteTempFiles();
		return res;
    }
    
//    public static double getFisherOneTailedGreater(int a, int b, int c, int d) {
//    	double p = 0.0;
//    	int z = Math.min(b, c);
//    	for(int i = 0; i <= z; i++) {
//    		p += getFisherP(a+i, b-i, c-i, d+i);
//    	}
//    	return p;
//    }
    
    public double getBinomialBackground(ContingencyTable ct, double bg_freq) {
    	
    	//add pseudo counts
    	int case_aff = ct.getA();
    	int n_cases = ct.getA()+ct.getB() +1;
    	int control_aff = ct.getC();
    	int n_controls = ct.getC() + ct.getD()+1;
    	int exp_cases = (int)(bg_freq * n_cases);
    	int exp_controls = (int)(bg_freq * n_controls);
    	
//    	NormalDistribution nd = new NormalDistribution();		
//		double p_val_case_vs_bg = 1 - new BinomialDistribution(n_cases+1, bg_freq).cumulativeProbability(case_aff);
//		double z_score_case_vs_bg = nd.inverseCumulativeProbability(p_val_case_vs_bg);
//		double p_val_control_vs_bg = 1 - new BinomialDistribution(n_controls+1, bg_freq).cumulativeProbability(control_aff);
//		double z_score_control_vs_bg = nd.inverseCumulativeProbability(p_val_control_vs_bg);
//		double z_score_diff = z_score_case_vs_bg - z_score_control_vs_bg;
//		return nd.cumulativeProbability(z_score_diff);
    	
    	code.clear();
    	if(exp_cases > case_aff) {
    		code.addRCode("z_case_bg<--qnorm(pbinom("+case_aff+","+n_cases+","+bg_freq+",lower.tail=T))");
    	} else {
    		code.addRCode("z_case_bg<-qnorm(pbinom("+case_aff+","+n_cases+","+bg_freq+",lower.tail=F))");
    	}
		
    	if(exp_controls > control_aff) {
    		code.addRCode("z_ctrl_bg<--qnorm(pbinom("+control_aff+","+n_controls+","+bg_freq+",lower.tail=T))");

    	} else {
    		code.addRCode("z_ctrl_bg<-qnorm(pbinom("+control_aff+","+n_controls+","+bg_freq+",lower.tail=F))");
    	}
		code.addRCode("res<-pnorm(z_case_bg-z_ctrl_bg)");

		r.runAndReturnResultOnline("res");
		double res = r.getParser().getAsDoubleArray("res")[0];
		r.deleteTempFiles();
		return res;
    }
    
    public double [] adjustP(double [] a, String method) {
    	code.clear();
    	
    	code.addDoubleArray("a", a);
    	code.addRCode("res<-p.adjust(a,\""+method+"\")");
    	r.runAndReturnResultOnline("res");
    	double [] res = r.getParser().getAsDoubleArray("res");
		r.deleteTempFiles();
    	return res;
    }
}

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
    
    public double getFisherOneTailedGreater(ContingencyTable ct) {
    	
    	//add pseudo counts
    	int x = ct.getA(); //+1;
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
    
//    public double getNormalApproximation(ContingencyTable ct, int pop_size, int pop_cases) {
//    	int N_case = ct.getA()+ct.getB() + pop_size;
//    	int D_case = ct.getA() + pop_cases; //has LOF
//    	int n_case = ct.getA()+ct.getB(); //has disease
//    	int x_case = ct.getA();
//    	double p_case = (double)D_case/(double)N_case;
//    	double z_case = (((double)x_case - 0.5) - (double)n_case * p_case)/(Math.sqrt((double)n_case*p_case*(1.0-p_case)));
//    	
//    	int N_ctrl = ct.getC() + ct.getD() + pop_size;
//    	int D_ctrl = ct.getC() + pop_cases;
//    	int n_ctrl = ct.getC() + ct.getD();
//    	int x_ctrl = ct.getC();
//    	double p_ctrl = (double)D_ctrl/(double)N_ctrl;
//    	double z_ctrl = (((double)x_ctrl - 0.5) - (double)n_ctrl * p_ctrl)/(Math.sqrt((double)n_ctrl*p_ctrl*(1.0-p_ctrl)));
//
//    	return z_case - z_ctrl;
//    }
    
    public double getBinomialBackground(ContingencyTable ct, double bg_freq) {
    	
    	//add pseudo counts
    	int case_aff = ct.getA();
    	int n_cases = ct.getA()+ct.getB() +1;
    	int control_aff = ct.getC();
    	int n_controls = ct.getC() + ct.getD()+1;
    	int exp_cases = (int)(bg_freq * n_cases);
    	int exp_controls = (int)(bg_freq * n_controls);
    	
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

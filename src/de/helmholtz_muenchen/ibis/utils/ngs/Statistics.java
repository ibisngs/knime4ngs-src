package de.helmholtz_muenchen.ibis.utils.ngs;

import rcaller.RCaller;
import rcaller.RCaller.FailurePolicy;
import rcaller.RCode;

public class Statistics {
	
	private RCaller r = null;
	private RCode code = null;
	
	public Statistics() {
		r = new RCaller();
		r.setRscriptExecutable("/home/software/bin/Rscript");
		r.setRExecutable("/home/software/bin/R");
		r.setFailurePolicy(FailurePolicy.CONTINUE);
		r.setMaxWaitTime(60000);
		code = new RCode();
		r.setRCode(code);
	}
	
	public void quit() {
		r.StopRCallerOnline();
		r.stopStreamConsumers();
	}
    
    public double [] getFisherOneTailedGreater(ContingencyTable [] ct) {
    	
    	int [] x = new int[ct.length];
    	int [] m = new int[ct.length];
    	int [] n = new int[ct.length];
    	int [] k = new int[ct.length];
    	
    	ContingencyTable my_table;
    	for(int i = 0; i < ct.length; i++) {
    		my_table = ct[i];
    		x[i] = my_table.getA(); //+1;
        	m[i] = x[i] + my_table.getB();
        	n[i] = my_table.getC() + my_table.getD();
        	k[i] = x[i] + my_table.getC();
    	}
    	
    	code.clear();
    	code.addIntArray("x", x);
    	code.addIntArray("m", m);
    	code.addIntArray("n", n);
    	code.addIntArray("k", k);
    	code.addRCode("res<-phyper(x-1,m,n,k,lower.tail=F)");
    	r.runAndReturnResultOnline("res");
    	double [] res = r.getParser().getAsDoubleArray("res");
		r.deleteTempFiles();
    	return res;
    }
    
public double [] getFisherTwoSided(ContingencyTable [] ct) {
    	
    	int [] a = new int[ct.length];
    	int [] b = new int[ct.length];
    	int [] c = new int[ct.length];
    	int [] d = new int[ct.length];
    	
    	ContingencyTable my_table;
    	for(int i = 0; i < ct.length; i++) {
    		my_table = ct[i];
    		a[i] = my_table.getA();
        	b[i] = my_table.getB();
        	c[i] = my_table.getC();
        	d[i] = my_table.getD();
    	}
    	
    	code.clear();
    	code.addIntArray("a", a);
    	code.addIntArray("b", b);
    	code.addIntArray("c", c);
    	code.addIntArray("d", d);
    	code.addRCode("df<- data.frame(a,b,c,d)");
    	code.addRCode("res<-apply(df,1,function(x) fisher.test(matrix(x,2))$p.value)");
    	r.runAndReturnResultOnline("res");
    	double [] res = r.getParser().getAsDoubleArray("res");
		r.deleteTempFiles();
    	return res;
    }
    
    public double [] getWilcoxon(ContingencyTable [] ct) {
    	
    	double [] result = new double[ct.length];
    	double [] res;
    	
    	for(int i = 0; i < ct.length; i++) {
    		ContingencyTable t = ct[i];
    		code.clear();
    		code.addRCode("case<-c(rep(1,"+t.getA()+"),rep(0,"+t.getB()+"))");
    		code.addRCode("ctrl<-c(rep(1,"+t.getC()+"),rep(0,"+t.getD()+"))");
    		code.addRCode("res<-wilcox.test(case,ctrl)$p.value");
    		r.runAndReturnResultOnline("res");
        	res = r.getParser().getAsDoubleArray("res");
    		r.deleteTempFiles();
        	result[i] = res[0];
    	}
    	
    	return result;
    }
    
    /*
     * returns same results as getHypergeometricBackground
     */
    public double [] getFisherExactBackground(ContingencyTable [] ct, int pop_size, double [] frequencies) {
    	int [] aff_case = new int[ct.length];
    	int [] aff_ctrl = new int[ct.length];
    	int [] un_case = new int[ct.length];
    	int [] un_ctrl = new int[ct.length];
    	int [] aff_pop = new int[ct.length];
    	int [] un_pop = new int[ct.length];
    	
    	ContingencyTable my_table;
    	for(int i = 0; i < ct.length; i++) {
    		my_table = ct[i];
    		aff_case[i] = my_table.getA()+1;
    		aff_ctrl[i] = my_table.getC()+1;	
    		un_case[i] = my_table.getB();
    		un_ctrl[i] = my_table.getD();
    		aff_pop[i] = (int)Math.round(frequencies[i]*pop_size);
    		un_pop[i] = pop_size - aff_pop[i];
    	}
    	
    	code.clear();
    	code.addIntArray("aff_case", aff_case);
    	code.addIntArray("aff_ctrl", aff_ctrl);
    	code.addIntArray("un_case", un_case);
    	code.addIntArray("un_ctrl", un_ctrl);
    	code.addIntArray("aff_pop", aff_pop);
    	code.addIntArray("un_pop", un_pop);
    	code.addRCode("df1 <- data.frame(aff_case,un_case,aff_pop,un_pop)");
    	code.addRCode("df2 <- data.frame(aff_ctrl,un_ctrl,aff_pop,un_pop)");
    	code.addRCode("z_case_bg<-qnorm(apply(df1,1,function(x) fisher.test(matrix(x,2),alternative=\"greater\")$p.value))");
    	code.addRCode("z_ctrl_bg<-qnorm(apply(df2,1,function(x) fisher.test(matrix(x,2),alternative=\"greater\")$p.value))");
    	code.addRCode("res<-pnorm(z_case_bg-z_ctrl_bg)");

		r.runAndReturnResultOnline("res");
		double [] res = r.getParser().getAsDoubleArray("res");
		r.deleteTempFiles();
		
		for(int i = 0 ; i< res.length; i++) {
			if(Double.compare(res[i], 0.0) == 0 || Double.compare(res[i], 1.0) == 0) {
				res[i] = Double.NaN;
			}
		}
		return res;
    }
    
    /*
     * returns the same results as getFisherExactBackground
     */
    public double [] getHypergeometricBackground(ContingencyTable [] ct, int pop_size, double [] frequencies) {
    	
    	int [] x_case = new int[ct.length];
    	int [] x_ctrl = new int[ct.length];
    	int [] m_case = new int[ct.length];
    	int [] m_ctrl = new int[ct.length];
    	int [] pop_cases = new int[ct.length];
    	
    	ContingencyTable my_table;
    	for(int i = 0; i < ct.length; i++) {
    		my_table = ct[i];
    		//add pseudo counts
    		x_case[i] = my_table.getA()+1;
    		x_ctrl[i] = my_table.getC()+1;	
    		m_case[i] = x_case[i] + my_table.getB();
    		m_ctrl[i] = x_ctrl[i] + my_table.getD();
    		pop_cases[i] = (int)Math.round(frequencies[i]*pop_size);
    	}
    	
//    	double pop_frac = (double)pop_cases/(double)pop_size;
//    	double case_frac = (double)x_case/(double)m_case;
//    	double ctrl_frac = (double)x_ctrl/(double)m_ctrl;
//    	if(pop_frac > case_frac) {
//        	code.addRCode("z_case_bg<--qnorm(phyper("+(x_case-1)+","+m_case+","+pop_size+","+(pop_cases+x_case)+",lower.tail=T))");
//    	} else {
//        	code.addRCode("z_case_bg<-qnorm(phyper("+(x_case-1)+","+m_case+","+pop_size+","+(pop_cases+x_case)+",lower.tail=F))");
//    	}
//    	
//    	if(pop_frac > ctrl_frac) {
//        	code.addRCode("z_ctrl_bg<--qnorm(phyper("+(x_ctrl-1)+","+m_ctrl+","+pop_size+","+(pop_cases+x_ctrl)+",lower.tail=T))");
//    	} else {
//        	code.addRCode("z_ctrl_bg<-qnorm(phyper("+(x_ctrl-1)+","+m_ctrl+","+pop_size+","+(pop_cases+x_ctrl)+",lower.tail=F))");
//    	}
    	
    	code.clear();
    	code.addIntArray("x_case", x_case);
    	code.addIntArray("x_ctrl", x_ctrl);
    	code.addIntArray("m_case", m_case);
    	code.addIntArray("m_ctrl", m_ctrl);
    	code.addIntArray("pop_cases", pop_cases);
    	code.addRCode("z_case_bg<-qnorm(phyper(x_case-1,m_case,"+pop_size+",pop_cases+x_case,lower.tail=F))");
    	code.addRCode("z_ctrl_bg<-qnorm(phyper(x_ctrl-1,m_ctrl,"+pop_size+",pop_cases+x_ctrl,lower.tail=F))");
    	code.addRCode("res<-pnorm(z_case_bg-z_ctrl_bg)");

		r.runAndReturnResultOnline("res");
		double [] res = r.getParser().getAsDoubleArray("res");
		r.deleteTempFiles();
		
		for(int i = 0 ; i< res.length; i++) {
			if(Double.compare(res[i], 0.0) == 0 || Double.compare(res[i], 1.0) == 0) {
				res[i] = Double.NaN;
			}
		}
		return res;
    }
    
    public double getWilcoxonSignedRankTest(ContingencyTable [] ct, int pop_size, double [] frequencies) {
    	
    	int [] x_case = new int[ct.length];
    	int [] x_ctrl = new int[ct.length];
    	int [] m_case = new int[ct.length];
    	int [] m_ctrl = new int[ct.length];
    	int [] pop_cases = new int[ct.length];
    	
    	ContingencyTable my_table;
    	for(int i = 0; i < ct.length; i++) {
    		my_table = ct[i];
    		//add pseudo counts
    		x_case[i] = my_table.getA();
    		x_ctrl[i] = my_table.getC();	
    		m_case[i] = x_case[i] + my_table.getB();
    		m_ctrl[i] = x_ctrl[i] + my_table.getD();
    		pop_cases[i] = (int)Math.round(frequencies[i]*pop_size);
    	}
    	
    	code.clear();
    	code.addIntArray("x_case", x_case);
    	code.addIntArray("x_ctrl", x_ctrl);
    	code.addIntArray("m_case", m_case);
    	code.addIntArray("m_ctrl", m_ctrl);
    	code.addIntArray("pop_cases", pop_cases);
//    	code.addRCode("case_bg<-phyper(x_case-1,m_case,"+pop_size+",pop_cases+x_case,lower.tail=F)");
//    	code.addRCode("ctrl_bg<-phyper(x_ctrl-1,m_ctrl,"+pop_size+",pop_cases+x_ctrl,lower.tail=F)");
    	code.addRCode("z_case_bg<-qnorm(phyper(x_case-1,m_case,"+pop_size+",pop_cases+x_case,lower.tail=F))");
    	code.addRCode("z_ctrl_bg<-qnorm(phyper(x_ctrl-1,m_ctrl,"+pop_size+",pop_cases+x_ctrl,lower.tail=F))");

    	code.addRCode("z_case_bg1<-replace(z_case_bg, !is.finite(z_case_bg) & z_case_bg==z_ctrl_bg,0)");
    	code.addRCode("z_ctrl_bg1<-replace(z_ctrl_bg, !is.finite(z_ctrl_bg) & z_case_bg==z_ctrl_bg,0)");
//    	code.addRCode("z_ctrl_bg1<-z_ctrl_bg[is.finite(z_case_bg) | is.finite(z_ctrl_bg)]");
    	
    	code.addRCode("res<-wilcox.test(z_case_bg1,z_ctrl_bg1,paired=TRUE)");
		code.addRCode("p<-res$p.value");
		r.setRCode(code);
		
		//debugging
//		if(set.equals("SECRETION_BY_CELL")) {
//			System.out.println("x_case");
//			System.out.println(Arrays.toString(x_case));
//			
//			r.runAndReturnResultOnline("z_case_bg");
//			System.out.println("z_case_bg");
//			System.out.println(Arrays.toString(r.getParser().getAsDoubleArray("z_case_bg")));
//			
//			System.out.println("x_ctrl");
//			System.out.println(Arrays.toString(x_ctrl));
//			r.runAndReturnResultOnline("z_ctrl_bg");
//			System.out.println("z_ctrl_bg");
//			System.out.println(Arrays.toString(r.getParser().getAsDoubleArray("z_ctrl_bg")));
//		}
		
		r.runAndReturnResultOnline("p");
		double res = r.getParser().getAsDoubleArray("p")[0];
		r.deleteTempFiles();
    	return res;

    }
    
    public double [] getBinomialBackground(ContingencyTable [] ct, double [] bg_freq, double pseudo_freq) {
    	
    	int [] case_aff = new int[ct.length];
    	int [] n_cases = new int[ct.length];
    	int [] control_aff = new int[ct.length];
    	int [] n_controls = new int[ct.length];
    	
    	ContingencyTable my_table;
    	for(int i = 0; i < ct.length; i++) {
    		my_table = ct[i];
    		//add pseudo counts
        	case_aff[i] = my_table.getA();
        	n_cases[i] = my_table.getA()+my_table.getB() +1;
        	control_aff[i] = my_table.getC();
        	n_controls[i] = my_table.getC() + my_table.getD()+1;
        	bg_freq[i] = (bg_freq[i] +pseudo_freq)/(1.0+pseudo_freq);
    	}
    	
//    	int exp_cases = (int)(bg_freq * n_cases);
//    	int exp_controls = (int)(bg_freq * n_controls);
//    	if(exp_cases > case_aff) {
//    		code.addRCode("z_case_bg<--qnorm(pbinom("+case_aff+","+n_cases+","+bg_freq+",lower.tail=T))");
//    	} else {
//    		code.addRCode("z_case_bg<-qnorm(pbinom("+case_aff+","+n_cases+","+bg_freq+",lower.tail=F))");
//    	}
//		
//    	if(exp_controls > control_aff) {
//    		code.addRCode("z_ctrl_bg<--qnorm(pbinom("+control_aff+","+n_controls+","+bg_freq+",lower.tail=T))");
//
//    	} else {
//    		code.addRCode("z_ctrl_bg<-qnorm(pbinom("+control_aff+","+n_controls+","+bg_freq+",lower.tail=F))");
//    	}
    	
    	code.clear();
    	code.addIntArray("case_aff", case_aff);
    	code.addIntArray("control_aff",control_aff);
    	code.addIntArray("n_cases", n_cases);
    	code.addIntArray("n_controls", n_controls);
    	code.addDoubleArray("bg_freq", bg_freq);
    	code.addRCode("z_case_bg<-qnorm(pbinom(case_aff,n_cases,bg_freq,lower.tail=F))");
    	code.addRCode("z_ctrl_bg<-qnorm(pbinom(control_aff,n_controls,bg_freq,lower.tail=F))");
		code.addRCode("res<-pnorm(z_case_bg-z_ctrl_bg)");

		r.runAndReturnResultOnline("res");
		double [] res = r.getParser().getAsDoubleArray("res");
		r.deleteTempFiles();
		
		for(int i = 0 ; i< res.length; i++) {
			if(Double.compare(res[i], 0.0) == 0 || Double.compare(res[i], 1.0) == 0) {
				res[i] = Double.NaN;
			}
		}
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

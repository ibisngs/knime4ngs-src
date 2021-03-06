/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.utils.ngs;

import java.util.List;
import java.util.Random;

import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.ngs.geneBasedAnalysis.CaseControlArray;
import rcaller.RCaller;
import rcaller.RCaller.FailurePolicy;
import rcaller.RCode;

public class Statistics {
	
	private RCaller r = null;
	private RCode code = null;
	private NodeLogger logger = null;
	
	public Statistics() {
		this(NodeLogger.getLogger(Statistics.class));
	}
	
	public Statistics(NodeLogger logger) {
		r = new RCaller();
		r.setRscriptExecutable("/home/software/bin/Rscript");
		r.setRExecutable("/home/software/bin/R");
//		r.setRscriptExecutable("/usr/bin/Rscript");
//		r.setRExecutable("/usr/bin/R");
		r.setFailurePolicy(FailurePolicy.CONTINUE);
		r.setMaxWaitTime(60000);
		code = new RCode();
		r.setRCode(code);
		this.logger = logger;
	}
	
	public void quit() {
		r.StopRCallerOnline();
		r.stopStreamConsumers();
	}
    
//    public double [] getFisherOneTailedGreater(ContingencyTable [] ct) {
//    	
//    	int [] x = new int[ct.length];
//    	int [] m = new int[ct.length];
//    	int [] n = new int[ct.length];
//    	int [] k = new int[ct.length];
//    	
//    	ContingencyTable my_table;
//    	for(int i = 0; i < ct.length; i++) {
//    		my_table = ct[i];
//    		x[i] = my_table.getA(); //+1;
//        	m[i] = x[i] + my_table.getB();
//        	n[i] = my_table.getC() + my_table.getD();
//        	k[i] = x[i] + my_table.getC();
//    	}
//    	
//    	code.clear();
//    	code.addIntArray("x", x);
//    	code.addIntArray("m", m);
//    	code.addIntArray("n", n);
//    	code.addIntArray("k", k);
//    	code.addRCode("res<-phyper(x-1,m,n,k,lower.tail=F)");
//    	r.runAndReturnResultOnline("res");
//    	double [] res = r.getParser().getAsDoubleArray("res");
//		r.deleteTempFiles();
//    	return res;
//    }
    
    public double [] getFisherTest(ContingencyTable [] ct, String alternative) {
    	
    	if(!(alternative.equals("greater") || alternative.equals("less"))) {
    		alternative = "two.sided";
    	}
    	
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
    	code.addRCode("res<-apply(df,1,function(x) fisher.test(matrix(x,2),alternative=\""+alternative+"\")$p.value)");
    	r.runAndReturnResultOnline("res");
    	double [] res = r.getParser().getAsDoubleArray("res");
		r.deleteTempFiles();
    	return res;
    }
    
    public double [] getWilcoxon(CaseControlArray [] cca, String alternative) {
    	
    	if(!(alternative.equals("greater") || alternative.equals("less"))) {
    		alternative = "two.sided";
    	}
    	
    	double [] result = new double[cca.length];
    	double [] res;
    	
//    	if(cca.length==0) return result;
//    	
//    	StringBuilder my_code = new StringBuilder("xaver<-list(list("+array2vector(cca[0].getCaseArray())+","+array2vector(cca[0].getControlArray())+")");
//    	for(CaseControlArray c: cca) {
//    		my_code.append(",list("+array2vector(c.getCaseArray())+","+array2vector(c.getControlArray())+")");
//    	}
//    	my_code.append(")");
//    	
//    	code.clear();
//    	code.addRCode(my_code.toString());
//    	code.addRCode("res<-sapply(xaver, function(x) wilcox.test(x[[1]],x[[2]])$p.value");
//    	r.runAndReturnResultOnline("res");
//    	res = r.getParser().getAsDoubleArray("res");
//		r.deleteTempFiles();
//    	return res;
    	
    	logger.debug(cca.length+" transcripts/genes to analyze");
    	for(int i = 0; i < cca.length; i++) {
    		if(i % 1000 == 0) {
    			logger.debug("Progress: "+i+"/"+cca.length);
    		}
    		CaseControlArray t = cca[i];
    		code.clear();
    		
    		code.addShortArray("case", t.getCaseArray());
    		code.addShortArray("control", t.getControlArray());
    		code.addRCode("res<-wilcox.test(case,control,alternative=\""+alternative+"\")$p.value");
    		r.runAndReturnResultOnline("res");
        	res = r.getParser().getAsDoubleArray("res");
    		r.deleteTempFiles();
        	result[i] = res[0];
    	}
    	
    	return result;
    }
    
//    private String array2vector(short [] array) {
//    	if(array.length==0) {
//    		return "c()";
//    	}
//    	StringBuilder res = new StringBuilder("("+array[0]);
//    	for(int i = 1; i < array.length; i++) {
//    		res.append(","+array[i]);
//    	}
//    	
//    	return res.toString();
//    }
    
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
    	code.addRCode("z_case_bg<-qnorm(phyper(x_case-1,m_case,"+pop_size+",pop_cases+x_case,lower.tail=F))");
    	code.addRCode("z_ctrl_bg<-qnorm(phyper(x_ctrl-1,m_ctrl,"+pop_size+",pop_cases+x_ctrl,lower.tail=F))");

    	code.addRCode("z_case_bg1<-replace(z_case_bg, !is.finite(z_case_bg) & z_case_bg==z_ctrl_bg,0)");
    	code.addRCode("z_ctrl_bg1<-replace(z_ctrl_bg, !is.finite(z_ctrl_bg) & z_case_bg==z_ctrl_bg,0)");
    	
    	code.addRCode("res<-wilcox.test(z_case_bg1,z_ctrl_bg1,paired=TRUE)");
		code.addRCode("p<-res$p.value");
		r.setRCode(code);
		r.runAndReturnResultOnline("p");
		double res = r.getParser().getAsDoubleArray("p")[0];
		r.deleteTempFiles();
    	return res;
    }
    
//    public double [] getBinomialBackground(ContingencyTable [] ct, double [] bg_freq, double pseudo_freq) {
//    	
//    	int [] case_aff = new int[ct.length];
//    	int [] n_cases = new int[ct.length];
//    	int [] control_aff = new int[ct.length];
//    	int [] n_controls = new int[ct.length];
//    	
//    	ContingencyTable my_table;
//    	for(int i = 0; i < ct.length; i++) {
//    		my_table = ct[i];
//    		//add pseudo counts
//        	case_aff[i] = my_table.getA();
//        	n_cases[i] = my_table.getA()+my_table.getB() +1;
//        	control_aff[i] = my_table.getC();
//        	n_controls[i] = my_table.getC() + my_table.getD()+1;
//        	bg_freq[i] = (bg_freq[i] +pseudo_freq)/(1.0+pseudo_freq);
//    	}
//    	
//    	code.clear();
//    	code.addIntArray("case_aff", case_aff);
//    	code.addIntArray("control_aff",control_aff);
//    	code.addIntArray("n_cases", n_cases);
//    	code.addIntArray("n_controls", n_controls);
//    	code.addDoubleArray("bg_freq", bg_freq);
//    	code.addRCode("z_case_bg<-qnorm(pbinom(case_aff,n_cases,bg_freq,lower.tail=F))");
//    	code.addRCode("z_ctrl_bg<-qnorm(pbinom(control_aff,n_controls,bg_freq,lower.tail=F))");
//		code.addRCode("res<-pnorm(z_case_bg-z_ctrl_bg)");
//
//		r.runAndReturnResultOnline("res");
//		double [] res = r.getParser().getAsDoubleArray("res");
//		r.deleteTempFiles();
//		
//		for(int i = 0 ; i< res.length; i++) {
//			if(Double.compare(res[i], 0.0) == 0 || Double.compare(res[i], 1.0) == 0) {
//				res[i] = Double.NaN;
//			}
//		}
//		return res;
//    }
    
    public double [] adjustP(double [] a, String method) {
    	code.clear();
    	
    	code.addDoubleArray("a", a);
    	code.addRCode("res<-p.adjust(a,\""+method+"\")");
    	r.runAndReturnResultOnline("res");
    	double [] res = r.getParser().getAsDoubleArray("res");
		r.deleteTempFiles();
    	return res;
    }
    
    /**
     * 
     * @param elements from which to choose
     * @param probs is a matrix of probabilities: each row corresponds to a vector of probability weights for obtaining the elements of the vector being sampled
     * @return
     */
    public int[][] getSamples(int [] elements, double[][]probs, int size) {
    	int [][] result = new int[probs.length][size];
    	code.clear();
    	code.addDoubleMatrix("probs", probs);
    	code.addIntArray("elements", elements);
    	code.addRCode("res<-apply(probs,1, function(x) sample(elements,"+size+", replace=TRUE, x))");
    	r.runAndReturnResultOnline("res");
    	
    	int[] res = r.getParser().getAsIntArray("res");
    	
    	for(int i = 0; i < res.length; i++) {
    		result[i/size][i%size] = res[i];
    	}
    	
    	r.deleteTempFiles();
		return result;
    }
    
//    public int[][] getSamples(int [] elements, List<double[]>probs_list, int size) {
//    	
//    	double[][] probs = new double[probs_list.size()][elements.length];
//    	for(int i = 0; i < probs_list.size(); i++) {
//    		probs[i] = probs_list.get(i);
//    	}
//    	
//    	return this.getSamples(elements, probs, size);
//    }
    
    public int[][] getSamples(int [] elements, List<double[]>probs_list, int size) {
    	
    	int [][] result = new int[probs_list.size()][size];
    	Random rand = new Random();
    	
    	double [] tmp_probs, intervals;
    	double tmp, sum;
    	for(int i = 0; i < probs_list.size(); i++) {
    		tmp_probs = probs_list.get(i);
    		
    		intervals = new double[tmp_probs.length];
    		sum = 0.0;
    		for(int k = 0; k < intervals.length; k++) {
    			intervals[k] = tmp_probs[k] + sum;
    			sum += tmp_probs[k];
    		}
    		
    		for(int j = 0; j < size; j++) {
    			tmp = rand.nextDouble();
    			for(int z = 0; z < intervals.length; z++) {
    				if(tmp < intervals[z]) {
    					result[i][j] = elements[z];
    					break;
    				}
    			}
    		}
    	}
    	
    	return result;
    }
    
//    public static void main(String [] args) {
//    	String [] gts = {"0/0","0/1","1/1"};
//    	double [][] probs = new double[][]{new double[]{0.8,0.15,0.05}, new double[]{0.05,0.15,0.8},new double[]{0.1,0.8,0.1}};
//    	Statistics stats = new Statistics();
//    	String [][] res = stats.getSamples(gts, probs,10);
//    	for(int i = 0; i < res.length; i++) {
//    		for(int j = 0; j < res[i].length; j++) {
//    			System.out.print(res[i][j]+" ");
//    		}
//    		System.out.println();
//    	}
//    	stats.quit();
//    }
}

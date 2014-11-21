package de.helmholtz_muenchen.ibis.ngs.depthofcoverage;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

public class DepthOfCoverage {

	
	public static void processCoverageFile(String DoCOutfile){
		
		int[] Coverage = new int[102];
		
	    try(BufferedReader br = new BufferedReader(new FileReader(DoCOutfile))) {
	       
	    	String line;
	    	double totalBasePairs = 0;
	    	
	    	br.readLine();	//Skip Header
	    	
	    	while ((line = br.readLine()) != null) {
	    		totalBasePairs++;
	    		int curr_coverage = Integer.parseInt(line.split("\t")[3]);

	    		if(curr_coverage>101){
	    			curr_coverage=101;
	    		}
	    		
	    		Coverage[curr_coverage]++;        	
	        }
	    	
		    //Everything counted, now get %-values
//		    double[] perc = new double[101];
//		    perc[0]=100;
	    	PrintWriter writer = new PrintWriter(DoCOutfile+"_DoCSummary.csv", "UTF-8");	    	
		    double remaining = totalBasePairs;
		    writer.print("100.00");
		    for(int i = 1; i<101;i++){
		    	remaining = remaining-Coverage[i-1];		//Anzahl von BasePairs mit Cov > Coverage[i]
		    	double value = remaining/totalBasePairs*100;
		    	value = (double)Math.round(value * 100) / 100;
//		    	perc[i] = value;
		    	writer.print("\t"+value);
		    }
		    writer.close();
		    
	    	
	    	
	    } catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
	    
	    
	}
	
}

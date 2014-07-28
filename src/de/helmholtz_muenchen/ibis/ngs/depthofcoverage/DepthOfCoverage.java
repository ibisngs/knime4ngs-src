package de.helmholtz_muenchen.ibis.ngs.depthofcoverage;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class DepthOfCoverage {

	
	public static void processCoverageFile(String DoCOutfile){
		
		int[] Coverage = new int[102];
		
	    try(BufferedReader br = new BufferedReader(new FileReader(DoCOutfile))) {
	       
	    	String line;
	    	int totalBasePairs = 0;
	    	while ((line = br.readLine()) != null) {
	    		totalBasePairs++;
	    		int curr_coverage = Integer.parseInt(line.split("\t")[3]);
	    		if(curr_coverage>101){
	    			curr_coverage=101;
	    		}
	    		Coverage[curr_coverage]++;        	
	        }

		    //Everything counted, now get %-values
		    double[] perc = new double[101];
		    int remaining = totalBasePairs;
		    for(int i = 0; i<101;i++){
		    	remaining = totalBasePairs-Coverage[i];		//Anzahl von BasePairs mit Cov > Coverage[i]
		    	perc[i] = remaining/totalBasePairs*100;
		    }
	    	
		    System.out.println(Arrays.toString(perc));
	    	
	    	
	    } catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
	    
	    
	}
	
}

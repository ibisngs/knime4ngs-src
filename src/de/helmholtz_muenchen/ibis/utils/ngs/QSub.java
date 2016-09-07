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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
@Deprecated
public class QSub {

	String memory = "";
	String prefixname = "";
	String q = "qsub -b y -j y -l vf=";//40G -N 'JMS Bowtie2' -o "
	
	public QSub(String com, String name, String memory, String logfile, Boolean wait) throws IOException, InterruptedException {
		
		String qs = q + memory + "G -N '" + name + "' -o '" + logfile + "' '" + com + "'";
		
    	System.out.println(qs);
    	
    	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", qs);
    	Process p = b.start();
    	p.waitFor();
    	
    	if(wait) {
	    	BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
	    	String s = null;
	    	StringBuffer nodeEntry = new StringBuffer(60);
	    	while ((s = stdInput.readLine()) != null) {
		    	nodeEntry.append(s+" ");
	    	}
	    	
	    	String[] flds = nodeEntry.toString().split(" ");    	
	    	String jobnumber = flds[2];
	    	
	    	System.out.println("Process " + jobnumber + " started.");
	    	
	    	while(true) {
	    		b = new ProcessBuilder("/bin/sh", "-c", "qstat");
	        	p = b.start();
	        	p.waitFor();
	        	
	        	stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
	        	s = null;
	        	nodeEntry = new StringBuffer(60);
	        	while ((s = stdInput.readLine()) != null) {
	    	    	nodeEntry.append(s+" ");
	        	}
	        	
	        	int idx = nodeEntry.toString().lastIndexOf(jobnumber);
	        	
	        	if(idx == -1) {
	        		break;
	        	}
	        	
	        	Thread.sleep(1000);
	        	
	    	}
	    	
	    	System.out.println("Process " + jobnumber + " finished.");
    	}
	}

	
}

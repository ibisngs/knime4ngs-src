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
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.SimpleDateFormat;
import java.util.Date;

public class ShowOutput {

	private static String logfile = "";
	
	
    public static StringBuffer getLogEntry(Process p, String call) throws IOException{
    	
    	String s = null;
    	BufferedReader stdInput = new BufferedReader(new 
                InputStreamReader(p.getInputStream()));

    	BufferedReader stdError = new BufferedReader(new 
                InputStreamReader(p.getErrorStream()));
      
    	StringBuffer nodeEntry = new StringBuffer(60);
    	nodeEntry.append(call+"\n");
    	nodeEntry.append("###Standard output of the command:###\n");
    	while ((s = stdInput.readLine()) != null) {
    	    	nodeEntry.append(s+"\n");
    	}
    	nodeEntry.append("###Standard error of the command:###\n");
    	while ((s = stdError.readLine()) != null) {
    		nodeEntry.append(s+"\n");
    	}
    	
    	return nodeEntry;
    }
    
    public static void writeLogFile(String entry){
    	try {
    	    BufferedWriter out = new BufferedWriter(new FileWriter(logfile, true));
     	    out.write(entry+"\n");
    	    out.close();
    	} catch (IOException e) { }
    }
    
    public static void writeLogFile(StringBuffer x){
    	try {
    	    BufferedWriter out = new BufferedWriter(new FileWriter(logfile, true));
     	    out.write(x.toString()+"\n");
    	    out.close();
    	} catch (IOException e) { }    	
    }
    
    public static String getNodeStartTime(String toolname){
    	String nodeStartTime = getCurrentTimeStamp();
    	return "Node "+toolname+" started at "+nodeStartTime+"\n";
    }
    public static String getNodeEndTime(){
    	String nodeEndTime = getCurrentTimeStamp();
    	String nodeEndString = "Node finished at "+nodeEndTime+"\n##########################################################################################\n##########################################################################################";
    	return nodeEndString;
    }
    
    public static String getCurrentTimeStamp() {
        SimpleDateFormat sdfDate = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");//dd/MM/yyyy
        Date now = new Date();
        String strDate = sdfDate.format(now);
        return strDate;
    }

    
    public static void setLogFile(String filepath){
    	if(logfile == "") {
    		logfile = filepath;
    	}
    }
    public static String getLogFile(){
    	return logfile;
    }
    
   
}


	
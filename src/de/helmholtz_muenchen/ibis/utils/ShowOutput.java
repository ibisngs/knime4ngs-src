package de.helmholtz_muenchen.ibis.utils;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class ShowOutput {

	private static String logfile = "";


	public static StringBuffer getLogEntryStdOut(Process p, String call) throws IOException{
		String s = null;
		StringBuffer nodeEntry = new StringBuffer(60);
		BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
		while ((s = stdInput.readLine()) != null) {
			nodeEntry.append(s+"\n");
		}
		return nodeEntry;
	}

	public static StringBuffer getLogEntryStdErr(Process p, String call) throws IOException{
		String s = null;
		StringBuffer nodeEntry = new StringBuffer(60);
		BufferedReader stdError = new BufferedReader(new InputStreamReader(p.getErrorStream()));
		while ((s = stdError.readLine()) != null) {
			nodeEntry.append(s+"\n");
		}
		return nodeEntry;
	}
	
	@Deprecated
	public static void writeLogFile(String entry){
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(logfile, true));
			out.write(entry+"\n");
			out.close();
		} catch (IOException e) { }
	}
	
	@Deprecated
	public static void writeLogFile(StringBuffer x){
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(logfile, true));
			out.write(x.toString()+"\n");
			out.close();
		} catch (IOException e) { }    	
	}

	public static String getNodeStartTime(String toolname){
		String nodeStartTime = Utils.getCurrentTimeStamp();
		return "Node "+toolname+" started at "+nodeStartTime+"\n";
	}
	public static String getNodeEndTime(){
		String nodeEndTime = Utils.getCurrentTimeStamp();
		String nodeEndString = "Node finished at "+nodeEndTime+"\n##########################################################################################\n##########################################################################################";
		return nodeEndString;
	}


	@Deprecated
	public static void setLogFile(String filepath){
		if(logfile == "") {
			logfile = filepath;
		}
	}
	@Deprecated
	public static String getLogFile(){
		return logfile;
	}


}



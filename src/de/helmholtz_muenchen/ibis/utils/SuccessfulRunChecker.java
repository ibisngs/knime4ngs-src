package de.helmholtz_muenchen.ibis.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.net.InetAddress;

public class SuccessfulRunChecker {
	
	private final BufferedWriter BW;
	
	public static final String LOCK_NAME = "lastRunCommand";
	public static final String LOCK_ENDING = ".klock";
	private static final String OK = "Terminated successfully!";
	
	/*public static void main(String args[]) throws IOException {
		
		boolean b = hasTerminatedSuccessfully(new File("/tmp/lock.test"), "test -r");
		System.out.println("termination state:" + b);
		
		if(!b) {
			SuccessfulRunChecker checker = new SuccessfulRunChecker(new File("/tmp/lock.test"), "test -r");
			checker.writeOK();
		}
	}*/
	
	
	public SuccessfulRunChecker(File lock, String command) throws IOException {
		if(lock != null) {
			if(!lock.getParentFile().exists()) 
				lock.getParentFile().mkdirs();
			
			BW = new BufferedWriter(new FileWriter(lock)); 
			writeCommand(command);
		}
		else
			BW = null;
	}
	
	private void writeCommand(String command) throws IOException {
		if(BW != null) {
			BW.write(command);
			BW.newLine();
			BW.write(InetAddress.getLocalHost().getHostName());
			BW.newLine();
			BW.flush();
		}
	}
	
	public void writeOK() throws IOException {
		if(BW != null) {
			BW.write(OK);
			BW.flush();
		}
	}
	
	public void finalize() {
		try  {
			if(BW != null)
				BW.close();
		}
		catch(Exception e) { e.printStackTrace(); }
	}
	
	public static boolean hasTerminatedSuccessfully(File lock, String command) {
		// check, if file exists
		if(!(lock != null && lock.isFile() && lock.exists()))
			return false;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(lock));
			String line;
			if((line = br.readLine()) != null) {
				// check, if command is the same
				if(line.equals(command)) {
					br.readLine();
					if((line = br.readLine()) != null) {
						if(OK.equals(line)) 
							return true;
					}
				}
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		return false;
	}
}

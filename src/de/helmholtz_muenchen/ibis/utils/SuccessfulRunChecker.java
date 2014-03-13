package de.helmholtz_muenchen.ibis.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class SuccessfulRunChecker {
	
	private final BufferedWriter BW;
	
	private static final String OK = "Terminated successfully!";
	
	/*public static void main(String args[]) throws IOException {
		boolean b = hasTerminatedSuccessfully(new File("/tmp/lock.test"), "test -r");
		System.out.println("terminateion sate:" + b);
		
		if(!b) {
			SuccessfulRunChecker checker = new SuccessfulRunChecker(new File("/tmp/lock.test"));
			checker.writeCommand("test -r");
			checker.writeOK();
		}
	}*/
	
	
	public SuccessfulRunChecker(File lock) throws IOException {
		BW = new BufferedWriter(new FileWriter(lock)); 
	}
	
	public void writeCommand(String command) throws IOException {
		BW.write(command);
		BW.newLine();
		BW.flush();
	}
	
	public void writeOK() throws IOException {
		BW.write(OK);
		BW.flush();
	}
	
	public void finalize() {
		try  {
			BW.close();
		}
		catch(Exception e) { e.printStackTrace(); }
	}
	
	public static boolean hasTerminatedSuccessfully(File lock, String command) {
		// check, if file exists
		if(!(lock.isFile() && lock.exists()))
			return false;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(lock));
			String line;
			if((line = br.readLine()) != null) {
				// check, if command is the same
				if(line.equals(command)) {
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

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
	
	/*examplary integration in execute method
	
	File lockFile = new File(...);
	String lockCommand = ...;
	boolean b = hasTerminatedSuccessfully(lockFile, lockCommand);
		
	if(!b) {
		SuccessfulRunChecker checker = new SuccessfulRunChecker(lockFile, lockCommand);
		int exitCode = Executor.executeCommandWithExitCode(...);
		if(exitCode==0) {
			checker.writeOK();
			checker.finalize();
		}
		
	}
	*/
	
	
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
						if(OK.equals(line)) {
							br.close();
							return true;
						}
							
					}
				}
			}
			br.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		return false;
	}
}

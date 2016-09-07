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

import java.io.File;

public class BinaryHandler {

	public final static String TOOL_PATH = "/home/software/bin";
	
	@Deprecated
	public static String checkToolAvailability(String BinaryName) {
		return checkToolAvailability(BinaryName, TOOL_PATH);
	}
	
	public static String checkToolAvailability(String BinaryName, String dir){
		
		FileSearch f = new FileSearch();
		f.searchDirectory(new File(dir), BinaryName);
		f.getResult();

		String ToolPath = null;
		int count = f.getResult().size();
		if(count ==0){
		    System.out.println("\nNo result found!");
		}else{
		    System.out.println("\nFound " + count + " result!\n");
		    for (String matched : f.getResult()){
		    	ToolPath = matched;
		    	System.out.println("Found : " + matched);
		    }
		}
		
		return ToolPath;
	}
	
	
	
	
	
	
}

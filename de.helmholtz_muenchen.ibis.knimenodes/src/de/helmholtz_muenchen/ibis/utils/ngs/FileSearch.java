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

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.knime.core.node.InvalidSettingsException;
 
public class FileSearch {
 
  private String fileNameToSearch;
  private List<String> result = new ArrayList<String>();
 
  public String getFileNameToSearch() {
	return fileNameToSearch;
  }
 
  public void setFileNameToSearch(String fileNameToSearch) {
	this.fileNameToSearch = fileNameToSearch;
  }
 
  public List<String> getResult() {
	return result;
  }
 
 
  public void searchDirectory(File directory, String fileNameToSearch) {
 
	setFileNameToSearch(fileNameToSearch);
 
	if (directory.isDirectory()) {
	    searchSuffix(directory);
	} else {
	    System.out.println(directory.getAbsoluteFile() + " is not a directory!");
	}
 
  }
 
  private void searchSuffix(File file) {
 
	if (file.isDirectory()) {
	  System.out.println("Searching directory ... " + file.getAbsoluteFile());
 
            //do you have permission to read this directory?	
	    if (file.canRead()) {
		for (File temp : file.listFiles()) {
		    if (temp.isDirectory()) {
			searchSuffix(temp);
		    } else {
			if (temp.getName().endsWith(getFileNameToSearch())) {			
			    result.add(temp.getAbsoluteFile().toString());
		    }
 
		}
	    }
 
	 } else {
		System.out.println(file.getAbsoluteFile() + "Permission Denied");
	 }
      }
 
  }
 
	/**
	 * Finds all Files in a given Directory that end with suffix
	 * @param dataDirectory
	 * @param regex
	 * @param Files2Merge
	 * @param logger 
	 * @throws InvalidSettingsException 
	 */
	public static void searchWithV(String dataDirectory, String regex, LinkedList<String> Files2Merge) throws InvalidSettingsException{
		File root = new File( dataDirectory );
		boolean foundFile = false;
	    File[] list = root.listFiles();
	        for ( File f : list ) {
	            if ( f.isDirectory() ) {
	            	searchWithV( f.getAbsolutePath() , regex, Files2Merge);
	            }else {
	            	String filePath = f.getAbsolutePath();
	            	if(filePath.matches(regex)){			//Add if file ends with specified regex
	            		Files2Merge.add("-V "+filePath);//, new String[]{filePath,"0",null,null,"-1","0"}); 
	            		foundFile = true;
	            	}
	            }
	        } 
	        
	        if(!foundFile){
	        	throw new InvalidSettingsException("Found no file that matches the regex: "+regex+" in folder: "+dataDirectory);
	        }
	}
	
	/**
	 * Finds all Files in a given Directory that end with suffix
	 * @param dataDirectory
	 * @param regex
	 * @param Files2Merge
	 * @param logger 
	 */
	public static void search(String dataDirectory, String regex, LinkedList<String> Files2Merge){
		File root = new File( dataDirectory );
	    File[] list = root.listFiles();
	        for ( File f : list ) {
	            if ( f.isDirectory() ) {
	            	search( f.getAbsolutePath() , regex, Files2Merge);
	            }else {
	            	String filePath = f.getAbsolutePath();
	            	if(filePath.endsWith(regex)){			//Add if file ends with specified regex
	            		Files2Merge.add(filePath);//, new String[]{filePath,"0",null,null,"-1","0"}); 
	            	}
	            }
	        } 
	}
  
  
}
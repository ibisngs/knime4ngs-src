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
package de.helmholtz_muenchen.ibis.utils.lofs;

import java.nio.file.Files;
import java.nio.file.Paths;

public class PathProcessor {
	
	/** processes a file path and returns the path without file extension, i.e /home/folder/file.txt -> /home/folder/file
	 * 
	 * @param path: path to a file
	 * @return path without file extension
	 */
	
	public static String getBase (String path){
		
		//get parent directory of file
        String directory = Paths.get(path).getParent().toString();
        //get file name without parent directory
        String filename=Paths.get(path).getFileName().toString();
        //split file name at .
        String[] splitn=filename.split("\\.");
        String basefilename=splitn[0];
        for(int i=1; i<splitn.length-1;i++){
        	basefilename=basefilename+"."+splitn[i];
        }
        //base name contains file name without extension and parent directory
        // concatenate parent directory and base name of the file
        String base=Paths.get(directory, basefilename).toString();
		
		return base;
	}
	
	/** processes a file path and returns the path without file extension, i.e /home/folder/file.txt -> txt
	 * 
	 * @param path to a file
	 * @return file extension
	 */
	
	// returns file extension of the path
	public static String getExt (String path){
		
		//get parent directory of file
        //get file name without parent directory
        String filename=Paths.get(path).getFileName().toString();
        //split file name at .
        String[] splitn=filename.split("\\.");
        String fileextension=splitn[splitn.length-1];
        
        return fileextension;
	}
	
	/** creates path to output file by concatenating filebase format and toolextension to a file path like this: filebase.toolextension.format
	 * it also checks if the path already exists and adds a number X to it: filebase.toolextensionX.format 
	 * 
	 * @param filebase output file without file extension
	 * @param format file extension of the output file
	 * @param toolextension String that is added to the file path which indicates the tool that creates the file
	 * @return 
	 * @throws Exception: prevents method from getting stuck in a while loop
	 */
	
	public static String createOutputFile(String filebase, String format, String toolextension) throws Exception{
		
		//path to output file
		String path =filebase+"."+toolextension+"."+format;
		
		//file does not exist
		if(!Files.exists(Paths.get(path))){
			return path;
		}
		
		//file exists
		else{
			
//			//add integer to file name and increment integer until the file does not exist
//			int n=1;
//			while(Files.exists(Paths.get(filebase+"."+toolextension+n+"."+format))){
//				n++;
//				
//				if(n==10000){
//					throw new Exception("Oops, I'm so sorry, something went wrong, there are too many files");
//				}
//			}
//			return filebase+"."+toolextension+n+"."+format;
			return path;
		}	
	}
	
	
}

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

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;

public class  Writer_Output {

	public String dir;
	public BufferedWriter w;
	
	
	/**
	 * create a buffered writer for a file 
	 * @param path: file to write
	 */
	
	public Writer_Output (String path){
		dir=path;
		initiateWriter(path);
	}
	
	public void initiateWriter (String path){
		try{
			Charset c = Charset.forName("UTF-8");
			w = Files.newBufferedWriter(Paths.get(path), c);
		}
		catch(IOException e){
			System.out.println("error creating the writer");
		}
	}
	
	/**
	 * writes string to file
	 * @param s: String which is written
	 */
	
	public void writeFile(String s){
		try{
			w.write(s);
		}
		catch(IOException e){
			System.out.println("error writing the file");
			
		}
	}
	
	/**
	 * writes String to file and adds a newline (like Sys.out.println)
	 * @param s: String which is written
	 */
	
	public void writeFileln(String s){
		writeFile(s+"\n");
	}
	
	/**
	 * closes the writer
	 */
	
	public void closew(){
		try{
			w.flush();
			w.close();
		}
		catch(IOException e){
			System.out.println("error closing the writer");
		}
	}
	
}

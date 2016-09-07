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
package de.helmholtz_muenchen.ibis.utils.threads;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class InputThread extends Thread{
	
    private InputStream input;		// InputStream generated from file -> read from STDIN
    private OutputStream output;	// OutputStream connected to STDIN of process
    
    private String filepath;

    
    public InputThread(OutputStream output, String file){
    	
    	try {
			this.input=new FileInputStream(new File(file));
		} 
    	catch (FileNotFoundException e) {
			e.printStackTrace();
		}
    	
    	this.output=output;
    	this.filepath=file;
    	
    }
    
    public void run() {
    	
    	System.out.println("starting stdinthread for: "+filepath);
    	
    	
    	
        try {
        	
            byte[] b = new byte[131072];
            int read = 0;
            
            // read= -1 -> end of file
            while (read > -1) {	
                read = input.read(b, 0, b.length);
                if (read > -1) {
                    output.write(b, 0, read);
                }
            }
            
            System.out.println("Closing InputThread for: "+filepath);
            
            input.close();
            output.close();
        } 
       catch (IOException e) {
            e.printStackTrace();
       } 

    }

}

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

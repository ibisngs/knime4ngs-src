package de.helmholtz_muenchen.ibis.utils.lofs;

import java.io.*;
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

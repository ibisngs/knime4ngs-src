package de.helmholtz_muenchen.ibis.utils.lofs;

import java.io.*;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;

public class FileInputReader {
	
	private BufferedReader reader;
	
	/**
	 * creates a buffered reader
	 * @param path: file to read
	 */
	public FileInputReader(String path) throws IOException{
		Charset c = Charset.forName("UTF-8");
		reader = Files.newBufferedReader(Paths.get(path),c);

	}
	
	/**
	 * reads one line from the file
	 * @return
	 */
	public String read(){

			String res;
			try {
				res = reader.readLine();
				return res;
			} catch (IOException e) {
				e.printStackTrace();
			}

		return "";
	}
	
	/**
	 * closes the writer
	 * @throws IOException
	 */
	public void closer() throws IOException{
			reader.close();
	}


}

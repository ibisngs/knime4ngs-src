/**
 * @author hastreiter
 */


package de.helmholtz_muenchen.ibis.utils.threads;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class StreamThread extends Thread {

	/**Input Stream**/
	InputStream stream;
	File OutFile;

	StreamThread(InputStream stream, String OutFilePath) {
		this.stream = stream;

		this.OutFile = new File(OutFilePath);
		// if file doesnt exists, then create it
		if (!OutFile.exists()) {
			try {
				System.out.println("Created new file: "+OutFilePath);
				OutFile.createNewFile();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	@Override
	public void run() {
		try {

			/**
			 * Uses byte array for handling binary data
			 */
			OutputStream outstream = new BufferedOutputStream(new FileOutputStream(OutFile.getAbsoluteFile())); 

			int read = 0;
			byte[] bytes = new byte[32768];

			while ((read = stream.read(bytes)) != -1){
				outstream.write(bytes,0,read);
			}

			System.out.println("Closing StreamThread for File: "+OutFile.getAbsoluteFile());
			outstream.close();

		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

}

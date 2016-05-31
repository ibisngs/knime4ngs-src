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
import java.util.Arrays;

public class StreamThread extends Thread {

	/**Input Stream**/
	InputStream stream;
	File OutFile;
	private boolean fillSB;
	private StringBuffer sb;
	private boolean fillFile;

	StreamThread(InputStream stream, String OutFilePath) {
		this.stream = stream;
		this.fillSB = false;
		this.sb = null;
		this.fillFile = true;
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
	
	public StreamThread(InputStream stream, String OutFilePath, StringBuffer sb) {
		this.stream = stream;
		
		this.sb = sb;

		if(OutFilePath != null)
			{
				this.fillFile = true;
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
		else
			{
				this.fillFile = false;
			}
		
		if(this.sb != null)
			{
				this.fillSB = true;
			}
		else
			this.fillSB = false;
		
	}

	@Override
	public void run() {
		try {

			/**
			 * Uses byte array for handling binary data
			 */
			OutputStream outstream = null;
			
			if(this.fillFile){
				try{
					outstream = new BufferedOutputStream(new FileOutputStream(this.OutFile.getAbsoluteFile())); 
				}catch(IOException e){
					System.out.println("Failed to open file: "+this.OutFile.getAbsoluteFile());
				}
			}

			
			int read = 0;
			byte[] bytes = new byte[32768];

			while ((read = this.stream.read(bytes)) != -1){
				if(this.fillSB)
					{
						byte[] printByte = Arrays.copyOfRange(bytes, 0, read);
						this.sb.append(new String(printByte));
//						System.out.print(new String(printByte));
					}
				if(this.fillFile)
					{
						outstream.write(bytes,0,read);
					}
				if(!this.fillFile && !this.fillSB)
					{
//						System.out.println("Discarding " + read + "bytes due to no file / Buffer to write to...");
					}
			}

			if(this.fillFile)
			{
				System.out.println("Closing StreamThread for File: "+this.OutFile.getAbsoluteFile());
				outstream.close();
			}

		} catch (IOException ioe) {
			ioe.printStackTrace();
		} catch (NullPointerException e){
			e.printStackTrace();
		}
	}

}

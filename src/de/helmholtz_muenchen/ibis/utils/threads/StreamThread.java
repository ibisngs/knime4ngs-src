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
			
			if(this.fillFile)
				outstream = new BufferedOutputStream(new FileOutputStream(this.OutFile.getAbsoluteFile())); 

			int read = 0;
			byte[] bytes = new byte[32768];

			while ((read = this.stream.read(bytes)) != -1){
				if(this.fillSB)
					{
						this.sb.append(new String(bytes));
					}
				if(this.fillFile)
					{
						outstream.write(bytes,0,read);
					}
				if(!this.fillFile && !this.fillSB)
					{
//						System.out.println("Discarding " + read + "bytes due to no file / Buffer to write to...");
					}
				try {
					Thread.sleep(1);
				} catch (InterruptedException e) {
					System.out.println("interrupted ");
					
					while ((read = this.stream.read(bytes)) != -1){
						if(this.fillSB)
							{
								this.sb.append(new String(bytes));
							}
						if(this.fillFile)
							{
								outstream.write(bytes,0,read);
							}
					}
					
					outstream.flush();
					outstream.close();
				}
			}

			if(this.fillFile)
			{
				System.out.println("Closing StreamThread for File: "+this.OutFile.getAbsoluteFile());
				outstream.close();
			}

		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

}

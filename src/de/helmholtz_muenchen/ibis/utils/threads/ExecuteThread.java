/**
 * @author hastreiter
 */


package de.helmholtz_muenchen.ibis.utils.threads;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.concurrent.Callable;

import org.apache.commons.lang.StringUtils;
import org.knime.core.node.NodeLogger;


public class ExecuteThread implements Callable<Boolean> {

	private final String[] command;					// Command to execute
	private Process p;						// Process for execution
	private final String stdOutFile;				// Process stdOut Filepath
	private final String stdErrFile;				// Process stdErr Filepath
	private final StringBuffer stdOutStr;
	private final StringBuffer stdErrStr;
	private StreamThread stdErrStream;
	private StreamThread stdOutStream;
	private final NodeLogger LOGGER;
	
	/**
	 * ExecuteThread that reads the Error and StdOutstreams
	 * @param threadname
	 * @param command
	 * @param stdOutFile
	 * @param stdErrFile
	 */
	public ExecuteThread(String command[], NodeLogger logger, String stdOutFile, String stdErrFile, StringBuffer stdOutStr, StringBuffer stdErrStr) {
		this.command = command;

		this.stdErrFile=stdErrFile;
		this.stdOutFile=stdOutFile;
		this.stdOutStr=stdOutStr;
		this.stdErrStr=stdErrStr;
		
		this.LOGGER = logger;
	}

	// This is the entry point for thread.
	@Override
	public Boolean call() throws Exception {
		LOGGER.info("Running thread : " + command[0]);

		//Start the process
		p = Runtime.getRuntime().exec(command);

//	  	ProcessBuilder b = new ProcessBuilder(command);
//		p = b.start();
		

		//If Output is written to stdout/stderr
		if(this.stdOutFile!=null){
			stdOutStream = new StreamThread(p.getInputStream(),stdOutFile);
			stdOutStream.start();

		}
		if(this.stdErrFile!=null){//Start both streams
			stdErrStream = new StreamThread(p.getErrorStream(),stdErrFile);
			stdErrStream.start();
		}

		// WAIT FOR PROCESS TO BE FINISHED
		p.waitFor();

		if(this.stdErrStr != null){
			this.stdErrStr.append(getLogEntryStdErr(p));
		}
		if(this.stdOutStr != null){
			this.stdOutStr.append(getLogEntryStdOut(p));
		}
		// INTERRUPT STREAMS
		if(this.stdOutFile!=null){
			stdOutStream.interrupt();
		}
		if(this.stdErrFile!=null){
			stdErrStream.interrupt();
		}
		LOGGER.info("thread will now exit");

		return new Boolean(p.exitValue()==0);

	}

	public String getSTDOUT(){
		if(stdOutStr!=null){
			return(this.stdOutStr.toString());
		}
		if(p!=null){
			try {
				return(getLogEntryStdOut(p).toString());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return("");
	}

	public String getSTDERR(){
		if(stdErrStr!=null){
			return(this.stdErrStr.toString());
		}
		if(p!=null){
			try {
				return(getLogEntryStdErr(p).toString());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return("");
	}


	/**
	 * Stops the execution thread and the called command
	 * @throws IOException 
	 */
	public void cancel() {
		LOGGER.info("Stopping process...");
		p.destroy();	
	}

	public String getCommand(){
		return(StringUtils.join(this.command, " "));
	}

	public static StringBuffer getLogEntryStdOut(Process p) throws IOException{
		String s = null;
		StringBuffer nodeEntry = new StringBuffer(60);
		BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
		while ((s = stdInput.readLine()) != null) {
			nodeEntry.append(s+"\n");
		}
		return nodeEntry;
	}

	public static StringBuffer getLogEntryStdErr(Process p) throws IOException{
		String s = null;
		StringBuffer nodeEntry = new StringBuffer(60);
		BufferedReader stdError = new BufferedReader(new InputStreamReader(p.getErrorStream()));
		while ((s = stdError.readLine()) != null) {
			nodeEntry.append(s+"\n");
		}
		return nodeEntry;
	}
}

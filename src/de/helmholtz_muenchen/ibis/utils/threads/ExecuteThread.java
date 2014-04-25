/**
 * @author hastreiter
 */


package de.helmholtz_muenchen.ibis.utils.threads;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.concurrent.Callable;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.IO;


public class ExecuteThread implements Callable<Boolean> {
	public static final int NUM_LINES_STDOUT_STDERR = 100;

	private final String[] command;					// Command to execute
	private Process p;						// Process for execution
	private final String stdOutFile;				// Process stdOut Filepath
	private final String stdErrFile;				// Process stdErr Filepath
	private final String stdInFile;					// stdin Filepath
	private final StringBuffer stdOutStr;
	private final StringBuffer stdErrStr;
	private StreamThread stdErrStream;
	private StreamThread stdOutStream;
	private InputThread stdInStream;
	private final NodeLogger LOGGER;
	private final String[] ENVIRONMENT;
	
	/**
	 * ExecuteThread starts process to execute command in a new thread. STDOUT and STDERR are redirected to files or StringBuffers
	 * @param command command to be executed
	 * @param logger NodeLogger to write messages to
	 * @param stdOutFile if not null STDOUT is redirected to this file
	 * @param stdErrFile if not null STDERR is redirected to this file
	 * @param stdOutStr if not null STDOUT is redirected to this StringBuffer
	 * @param stdErrStr if not null STDERR is redirected to this StringBuffer
	 * @param enableEscape enables parameter escaping
	 */
	public ExecuteThread(String[] command, NodeLogger logger, String stdOutFile, String stdErrFile, String stdInFile, StringBuffer stdOutStr, StringBuffer stdErrStr, String[] Environment) {
		this.command = command;

		this.stdErrFile=stdErrFile;
		this.stdOutFile=stdOutFile;
		this.stdInFile=stdInFile;
		this.stdOutStr=stdOutStr;
		this.stdErrStr=stdErrStr;

		this.LOGGER = logger;
		
		// workaround for letting jar files run properly!
		if(Environment == null && command[0].startsWith("java -jar"))
			this.ENVIRONMENT = new String[0];
		else
			this.ENVIRONMENT = Environment;
	}


	/**
	 * Starts to run the command in a new thread and catches STDOUT and STDERR
	 */
	@Override
	public Boolean call() throws Exception {
		LOGGER.info("Running command: " + getCommand(this.command));
		
		//Start the process
		if(this.command.length==1){
			p = Runtime.getRuntime().exec(this.command[0], this.ENVIRONMENT);
		}else{
			p = Runtime.getRuntime().exec(this.command, this.ENVIRONMENT);
		}
		

		//If Output is written to stdout/stderr
		if(this.stdOutFile!=null){
			stdOutStream = new StreamThread(p.getInputStream(),stdOutFile);
			stdOutStream.start();

		}
		if(this.stdErrFile!=null){//Start both streams
			stdErrStream = new StreamThread(p.getErrorStream(),stdErrFile);
			stdErrStream.start();
		}
		if(this.stdInFile!=null){
			stdInStream = new InputThread(p.getOutputStream(), stdInFile);
			stdInStream.start();
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
		if(this.stdInFile!=null){
			stdInStream.interrupt();
		}
		LOGGER.info("finished command "+ command[0]);

		return new Boolean(p.exitValue()==0);
	}

	/**
	 * Returns the last x lines from STDOUT
	 * @return String containing STDOUT of the executed command
	 */
	public String getSTDOUT(){
		if(stdOutStr!=null){
			return(this.stdOutStr.toString());
		}
		if(this.stdOutFile != null){
			try {
				return(IO.tail(new File(this.stdOutFile), NUM_LINES_STDOUT_STDERR));
			} catch (IOException e) {
				this.LOGGER.warn("Can't read STDOUT file " + this.stdOutFile);
				e.printStackTrace();
			}
		}
		if(p!=null){
			try {
				return(getLogEntryStdOut(p).toString());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return("No STDOUT available!");
	}

	/**
	 * Returns the last x lines from STDERR
	 * @return String containing STDERR of the executed command
	 */
	public String getSTDERR(){
		if(this.stdErrStr!=null){
			return(this.stdErrStr.toString());
		}
		if(this.stdErrFile != null){
			try {
				return(IO.tail(new File(this.stdErrFile), NUM_LINES_STDOUT_STDERR));
			} catch (IOException e) {
				this.LOGGER.warn("Can't read STDERR file " + this.stdErrFile);
				e.printStackTrace();
			}
		}
		if(p!=null){
			try {
				return(getLogEntryStdErr(p).toString());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return("No STDERR available!");
	}
	
	public String printLogs(){
		return("LOGS OF "+ this.command[0] +"\nSTDOUT: \n" + this.getSTDOUT() + "\n\nSTDERR:\n" + this.getSTDERR());
	}
	/**
	 * Stops the execution thread and the called command
	 */
	public void cancel() {
		LOGGER.info("Stopping process...");
		p.destroy();	
	}

	/**
	 * Catch STDOUT of process and return it as StringBuffer
	 * @param p the process
	 * @return StringBuffer containing the STDOUT of the given process
	 * @throws IOException
	 */
	private static StringBuffer getLogEntryStdOut(Process p) throws IOException{
		String s = null;
		StringBuffer nodeEntry = new StringBuffer(60);
		BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
		while ((s = stdInput.readLine()) != null) {
			nodeEntry.append(s+"\n");
		}
		return nodeEntry;
	}

	/**
	 * Catch STDERR of process and return it as StringBuffer
	 * @param p the process
	 * @return StringBuffer containing the STDERR of the given process
	 * @throws IOException
	 */
	private static StringBuffer getLogEntryStdErr(Process p) throws IOException{
		String s = null;
		StringBuffer nodeEntry = new StringBuffer(60);
		BufferedReader stdError = new BufferedReader(new InputStreamReader(p.getErrorStream()));
		while ((s = stdError.readLine()) != null) {
			nodeEntry.append(s+"\n");
		}
		return nodeEntry;
	}
		
	public static String getCommand(String[] command){
		return StringUtils.join(command, "");
	}
	
	/**
	 * Gets the exit code
	 * @return
	 * @throws IllegalThreadStateException
	 */
	public int getExitCode() throws IllegalThreadStateException {
		return p.exitValue();
	}
}

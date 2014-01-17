/**
 * @author hastreiter
 */


package de.helmholtz_muenchen.ibis.utils.threads;

import java.io.IOException;
import java.lang.reflect.Field;


public class ExecuteThread implements Runnable {


	private String name; 					// name of thread
	private boolean finished; 				// Thread Status
	private Thread executingThread;			// The executing thread
	private String command;					// Command to execute
	private Process p;						// Process for execution
	private String stdOutFile;				// Process stdOut Filepath
	private String stdErrFile;				// Process stdErr Filepath


	/**
	 * ExecuteThread that reads the Error and StdOutstreams
	 * @param threadname
	 * @param command
	 * @param stdOutFile
	 * @param stdErrFile
	 */
	public ExecuteThread(String threadname, String command, String stdOutFile, String stdErrFile) {
		this.name = threadname;
		this.executingThread = new Thread(this, name);
		System.out.println("New execute thread: " + executingThread);
		this.finished = false;
		this.command = command;

		this.stdErrFile=stdErrFile;
		this.stdOutFile=stdOutFile;
	}


	// This is the entry point for thread.
	public void run() {
		try {
			System.out.println("Running command: "+command);

			//Start the process
			p = Runtime.getRuntime().exec(command);

			//If Output is written to stdout/stderr
			StreamThread stdErrStream = null;
			StreamThread stdOutStream = null;

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
			finished = true;

			// INTERRUPT STREAMS
			if(this.stdOutFile!=null){
				stdOutStream.interrupt();
			}
			if(this.stdErrFile!=null){
				stdErrStream.interrupt();
			}


		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println(this.name+" thread will now exit");
	}

	/**
	 * Stops the execution thread and the called command
	 */
	public void stop() {
		System.out.println("Stopping process...");
		p.destroy();
		//Close the Streams
		try {
			p.getErrorStream().close();
			p.getOutputStream().close();
			p.getInputStream().close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		try {
			p.waitFor();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}			    	
	}

	/**
	 * Starts the execution thread
	 */
	public void start(){
		executingThread.start();
	}

	public String getName(){
		return executingThread.getName();
	}

	/**
	 * Returns true if process finished
	 * @return
	 */
	public boolean isDone(){
		return this.finished;
	}

	public String getCommand(){
		return this.command;
	}

	/**
	 * Forces the main thread to wait until sub-threads finish
	 */
	public void waitUntilFinished(){
		try {
			executingThread.join();
		} catch (InterruptedException e) {
		}
	}

	@SuppressWarnings("unused")
	@Deprecated
	private void getPID(Process p){
		try {
			@SuppressWarnings("rawtypes")
			Class clazz = Class.forName("java.lang.UNIXProcess");
			Field pidField = clazz.getDeclaredField("pid");
			pidField.setAccessible(true);
			Object value = pidField.get(p);
			System.err.println("pid = " + value);
		} catch (Throwable e) {
			e.printStackTrace();
		}
	}
}

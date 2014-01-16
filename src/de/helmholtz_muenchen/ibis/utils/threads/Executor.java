/**
 * @author hastreiter
 */

package de.helmholtz_muenchen.ibis.utils.threads;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.knime.core.node.ExecutionContext;

public class Executor {
	
	private boolean run;
	
	/**
	 * Executor object that handles the tool executions.  
	 */
	public Executor(){
		this.run = true;
	}

	/**
	 * Multithreaded execution of command line calls. Uses a monitoring thread for checking if user canceled the execution thread.
	 * @param command : The command string. Note: No /bin/sh is used, therefore do not use ">" for directing output into files, use the tool flag or executeCommand with StdOut/StdError support instead 
	 * @param CommandName : String that will be used to name the threads 
	 * @param exec : The Knime ExecutionContext
	 */
	public void executeCommand(String command,String CommandName, ExecutionContext exec){
	 if(run){
		//Create Threadpool + the execution and monitoring threads 
		ExecutorService pool = Executors.newFixedThreadPool(2);
		ExecuteThread ex = new ExecuteThread(CommandName+" Executor",command);
		MonitorThread m = new MonitorThread(CommandName+" Monitor", ex,exec);

		//Starts both threads and waits for successful termination
		pool.execute(ex);
		pool.execute(m);	
	
		//Now wait for them
		ex.waitUntilFinished();
		m.waitUntilFinished();
		pool.shutdown();
        while (!pool.isTerminated()) {
        }
		System.out.println("ThreadPool closed.");
		
		//Check if execution was canceled to avoid follow-up processes
		if(m.isExecutionCanceled()){
			this.run = false;
		}
		
	 }else{
		System.out.println("Skipping command: "+command+"\nExcution was canceled..");
	 }
	}
	
	
	
	
	
	/**
	 * Multithreaded execution of command line calls. Uses a monitoring thread for checking if user canceled the execution thread.
	 * @param command : The command string. Note: No /bin/sh is used, therefore do not use ">" for directing output into files. Output to Stdout will be captured
	 * @param CommandName : String that will be used to name the threads 
	 * @param exec : The Knime ExecutionContext
	 * @param stdOutFile : File in which stdout will be written
	 */
		public void executeCommand(String command,String CommandName, ExecutionContext exec, String stdOutFile){
			 if(run){
				//Create Threadpool + the execution and monitoring threads 
				ExecutorService pool = Executors.newFixedThreadPool(2);
				ExecuteThread ex = new ExecuteThread(CommandName+" Executor",command, stdOutFile);
				MonitorThread m = new MonitorThread(CommandName+" Monitor", ex,exec);

				//Starts both threads and waits for successful termination
				pool.execute(ex);
				pool.execute(m);	
			
				//Now wait for them
				ex.waitUntilFinished();
				m.waitUntilFinished();
				pool.shutdown();
		        while (!pool.isTerminated()) {
		        }
				System.out.println("ThreadPool closed.");
				
				//Check if execution was canceled to avoid follow-up processes
				if(m.isExecutionCanceled()){
					this.run = false;
				}
				
			 }else{
				System.out.println("Skipping command: "+command+"\nExcution was canceled..");
			 }
			}
	
	
	
	
	
	
	
/**
 * Multithreaded execution of command line calls. Uses a monitoring thread for checking if user canceled the execution thread.
 * @param command : The command string. Note: No /bin/sh is used, therefore do not use ">" for directing output into files. Output to Stdout and Stderr will be captured
 * @param CommandName : String that will be used to name the threads 
 * @param exec : The Knime ExecutionContext
 * @param stdOutFile : File in which stdout will be written
 * @param stdErrFile : File in which stderr will be written
 */
	public void executeCommand(String command,String CommandName, ExecutionContext exec, String stdOutFile, String stdErrFile){
		 if(run){
			//Create Threadpool + the execution and monitoring threads 
			ExecutorService pool = Executors.newFixedThreadPool(2);
			ExecuteThread ex = new ExecuteThread(CommandName+" Executor",command, stdOutFile, stdErrFile);
			MonitorThread m = new MonitorThread(CommandName+" Monitor", ex,exec);

			//Starts both threads and waits for successful termination
			pool.execute(ex);
			pool.execute(m);	
		
			//Now wait for them
			ex.waitUntilFinished();
			m.waitUntilFinished();
			pool.shutdown();
	        while (!pool.isTerminated()) {
	        }
			System.out.println("ThreadPool closed.");
			
			//Check if execution was canceled to avoid follow-up processes
			if(m.isExecutionCanceled()){
				this.run = false;
			}
			
		 }else{
			System.out.println("Skipping command: "+command+"\nExcution was canceled..");
		 }
		}
	
	
}

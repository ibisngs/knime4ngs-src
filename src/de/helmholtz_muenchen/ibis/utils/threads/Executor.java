/**
 * @author hastreiter
 */

package de.helmholtz_muenchen.ibis.utils.threads;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.NodeLogger;

public class Executor {


	/**
	 * Multithreaded execution of command line calls. Uses a monitoring thread for checking if user canceled the execution thread.
	 * @param command : The command string. Note: No /bin/sh is used, therefore do not use ">" for directing output into files, use the tool flag or executeCommand with StdOut/StdError support instead 
	 * @param CommandName : String that will be used to name the threads 
	 * @param exec : The Knime ExecutionContext
	 * @throws Exception 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, NodeLogger logger) throws Exception{
		Executor.executeCommand(command, exec, logger, null, null, null, null);
	}

	/**
	 * Multithreaded execution of command line calls. Uses a monitoring thread for checking if user canceled the execution thread.
	 * @param command : The command string. Note: No /bin/sh is used, therefore do not use ">" for directing output into files. Output to Stdout will be captured
	 * @param CommandName : String that will be used to name the threads 
	 * @param exec : The Knime ExecutionContext
	 * @param stdOutFile : File in which stdout will be written
	 * @throws Exception 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, NodeLogger logger, String stdOutFile) throws Exception{
		Executor.executeCommand(command, exec, logger, stdOutFile, null, null, null);
	}
	
	public static void executeCommand(String[] command, ExecutionContext exec, NodeLogger logger, StringBuffer stdOut, StringBuffer stdErr) throws Exception{
		Executor.executeCommand(command, exec, logger, null, null, stdOut, stdErr);
	}


	/**
	 * Multithreaded execution of command line calls. Uses a monitoring thread for checking if user canceled the execution thread.
	 * @param command : The command string. Note: No /bin/sh is used, therefore do not use ">" for directing output into files. Output to Stdout and Stderr will be captured
	 * @param CommandName : String that will be used to name the threads 
	 * @param exec : The Knime ExecutionContext
	 * @param stdOutFile : File in which stdout will be written
	 * @param stdErrFile : File in which stderr will be written
	 * @throws CanceledExecutionException 
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, NodeLogger logger, String stdOutFile, String stdErrFile, StringBuffer stdOut, StringBuffer stdErr) throws Exception {
			//Create Threadpool + the execution and monitoring threads 
			ExecutorService pool = Executors.newSingleThreadExecutor();
			
			ExecuteThread executorThread   = new ExecuteThread(command, logger, stdOutFile, stdErrFile, stdOut, stdErr);
			Future<Boolean> executorResult = pool.submit(executorThread);
			
			// wait for executorThread to finish and meanwhile monitor if node was cancled
			while (!executorResult.isDone()) {
				// if cancel was requested, an exception will be thrown
				try{
					exec.checkCanceled();
				} catch (CanceledExecutionException e){
					// kill jobs
					executorThread.cancel();
					pool.shutdownNow();
					while (!pool.isTerminated()) {
					}
					throw(new CanceledExecutionException("Canceled Execution of command " + command[0] + "\nSTDERR:\n" + executorThread.getSTDERR() + "\n\nSTDOUT:\n" +  executorThread.getSTDOUT()));
				}
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {
					// ignore interrupted exception
				}
			}

			// check if finished successfully
			Boolean finishedSuccessfully = executorResult.get();
			if(!finishedSuccessfully){
				throw(new CanceledExecutionException("Could not successfully finish command " + command[0] + "\nSTDERR:\n" + executorThread.getSTDERR() + "\n\nSTDOUT:\n" +  executorThread.getSTDOUT()));
			}
			logger.info("ThreadPool closed.");
	}

}

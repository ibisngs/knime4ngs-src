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
	 * Execute a command in new thread, redirect its STDOUT and STDERR and 
	 * cancel the command if it is requested by the user via the given ExecutionContext
	 * STDOUT and STDERR are omitted in this case
	 * @param command String array contains the program/path to the script in the first entry and all commandline parameters in the following entries
	 * @param exec ExecutionContext
	 * @param logger a NodeLogger to write messages to
	 * @throws CanceledExecutionException 
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, NodeLogger logger) throws CanceledExecutionException, InterruptedException, ExecutionException {
		Executor.executeCommand(command, exec, null, logger, null, null, null, null);
	}

	/**
	 * Execute a command in new thread, redirect its STDOUT and STDERR and 
	 * cancel the command if it is requested by the user via the given ExecutionContext
	 * STDERR is omitted in this case
	 * @param command String array contains the program/path to the script in the first entry and all commandline parameters in the following entries
	 * @param exec ExecutionContext
	 * @param logger a NodeLogger to write messages to
	 * @param stdOutFile path to file to write STDOUT to, omitted if null
	 * @throws CanceledExecutionException 
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, NodeLogger logger, String stdOutFile) throws CanceledExecutionException, InterruptedException, ExecutionException{
		Executor.executeCommand(command, exec, null, logger, stdOutFile, null, null, null);
	}


	/**
	 * Execute a command in new thread, redirect its STDOUT and STDERR and 
	 * cancel the command if it is requested by the user via the given ExecutionContext
	 * STDOUT and STDERR are written to StringBuffers
	 * @param command String array contains the program/path to the script in the first entry and all commandline parameters in the following entries
	 * @param exec ExecutionContext
	 * @param logger a NodeLogger to write messages to
	 * @param stdOut  StringBuffer to write STDOUT to, omitted if null
	 * @param stdErr  StringBuffer to write STDERR to, omitted if null
	 * @throws CanceledExecutionException 
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, NodeLogger logger, StringBuffer stdOut, StringBuffer stdErr) throws CanceledExecutionException, InterruptedException, ExecutionException{
		Executor.executeCommand(command, exec, null, logger, null, null, stdOut, stdErr);
	}
	
	/**
	 * Execute a command in new thread, redirect its STDOUT and STDERR and 
	 * cancel the command if it is requested by the user via the given ExecutionContext
	 * STDOUT and STDERR are written to StringBuffers
	 * @param command String array contains the program/path to the script in the first entry and all commandline parameters in the following entries
	 * @param exec ExecutionContext
	 * @param logger a NodeLogger to write messages to
	 * @param stdOut  StringBuffer to write STDOUT to, omitted if null
	 * @param stdErr  StringBuffer to write STDERR to, omitted if null
	 * @throws CanceledExecutionException 
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, String[] environment, NodeLogger logger, StringBuffer stdOut, StringBuffer stdErr) throws CanceledExecutionException, InterruptedException, ExecutionException{
		Executor.executeCommand(command, exec, environment, logger, null, null, stdOut, stdErr);
	}


	/**
	 * Execute a command in new thread, redirect its STDOUT and STDERR and 
	 * cancel the command if it is requested by the user via the given ExecutionContext
	 * @param command String array contains the program/path to the script in the first entry and all commandline parameters in the following entries
	 * @param exec ExecutionContext
	 * @param logger a NodeLogger to write messages to
	 * @param stdOutFile path to file to write STDOUT to, omitted if null
	 * @param stdErrFile path to file to write STDERR to, omitted if null
	 * @param stdOut  StringBuffer to write STDOUT to, omitted if null
	 * @param stdErr  StringBuffer to write STDERR to, omitted if null
	 * @throws CanceledExecutionException 
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 * @throws Exception
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, String[] environment, NodeLogger logger, String stdOutFile, String stdErrFile, StringBuffer stdOut, StringBuffer stdErr) throws CanceledExecutionException, InterruptedException, ExecutionException  {
		exec.checkCanceled();
		//Create Threadpool + the execution and monitoring threads 
		ExecutorService pool = Executors.newSingleThreadExecutor();

		ExecuteThread executorThread   = new ExecuteThread(command, logger, stdOutFile, stdErrFile, stdOut, stdErr, environment);
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
				throw(new CanceledExecutionException("Canceled Execution of command " + ExecuteThread.getCommand(command)));
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
			logger.error(executorThread.printLogs());
			
			throw(new CanceledExecutionException("Could not successfully finish command " + ExecuteThread.getCommand(command)));
		}
		//logger.info("ThreadPool closed successfully.");
	}
}

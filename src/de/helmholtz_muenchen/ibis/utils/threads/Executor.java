/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package de.helmholtz_muenchen.ibis.utils.threads;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.IO;

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
	 * @throws UnsuccessfulExecutionException 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, NodeLogger logger) throws CanceledExecutionException, InterruptedException, ExecutionException, UnsuccessfulExecutionException {
		Executor.executeCommand(command, exec, null, logger, null, null, null, null, null,null);
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
	 * @throws UnsuccessfulExecutionException 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, NodeLogger logger, String stdOutFile) throws CanceledExecutionException, InterruptedException, ExecutionException, UnsuccessfulExecutionException{
		Executor.executeCommand(command, exec, null, logger, stdOutFile, null, null, null, null,null);
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
	 * @throws InterruptedException
	 * @throws ExecutionException
	 * @throws UnsuccessfulExecutionException
	 */
	public static void executeCommand(String[] command, ExecutionContext exec,String[] environment, NodeLogger logger, String stdOutFile) throws CanceledExecutionException, InterruptedException, ExecutionException, UnsuccessfulExecutionException{
		Executor.executeCommand(command, exec, environment, logger, stdOutFile, null, null, null, null,null);
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
	 * @throws UnsuccessfulExecutionException 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, NodeLogger logger, StringBuffer stdOut, StringBuffer stdErr) throws CanceledExecutionException, InterruptedException, ExecutionException, UnsuccessfulExecutionException{
		Executor.executeCommand(command, exec, null, logger, null, null, stdOut, stdErr, null,null);
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
	 * @throws UnsuccessfulExecutionException 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, String[] environment, NodeLogger logger, StringBuffer stdOut, StringBuffer stdErr) throws CanceledExecutionException, InterruptedException, ExecutionException, UnsuccessfulExecutionException{
		Executor.executeCommand(command, exec, environment, logger, null, null, stdOut, stdErr, null,null);
	}
	
	/**
	 * Execute a command in new thread, redirect its STDOUT and STDERR and 
	 * cancel the command if it is requested by the user via the given ExecutionContext
	 * STDOUT and STDERR are written to files
	 * STDIN is read from file
	 * @param command String array contains the program/path to the script in the first entry and all commandline parameters in the following entries
	 * @param exec ExecutionContext
	 * @param logger a NodeLogger to write messages to
	 * @param stdOutFile  File to write STDOUT to, omitted if null
	 * @param stdErrFile  File to write STDERR to, omitted if null
	 * @param stdOut  StringBuffer to write STDOUT to, omitted if null
	 * @param stdErr  StringBuffer to write STDERR to, omitted if null
	 * @throws CanceledExecutionException 
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 * @throws UnsuccessfulExecutionException 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, String [] environment, NodeLogger logger, String stdOutFile, String stdErrFile, StringBuffer stdOut, StringBuffer stdErr) throws CanceledExecutionException, InterruptedException, ExecutionException, UnsuccessfulExecutionException{
		Executor.executeCommand(command, exec, environment, logger, stdOutFile, stdErrFile, stdOut, stdErr, null,null);
	}
	
	/**
	 * Execute a command in new thread, redirect its STDOUT and STDERR and 
	 * cancel the command if it is requested by the user via the given ExecutionContext
	 * STDOUT and STDERR are written to files
	 * STDIN is read from file
	 * @param command String array contains the program/path to the script in the first entry and all commandline parameters in the following entries
	 * @param exec ExecutionContext
	 * @param logger a NodeLogger to write messages to
	 * @param stdOutFile  File to write STDOUT to, omitted if null
	 * @param stdErrFile  File to write STDERR to, omitted if null
	 * @throws CanceledExecutionException 
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 * @throws UnsuccessfulExecutionException 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, String [] environment, NodeLogger logger, String stdOutFile, String stdErrFile) throws CanceledExecutionException, InterruptedException, ExecutionException, UnsuccessfulExecutionException{
		Executor.executeCommand(command, exec, environment, logger, stdOutFile, stdErrFile, null, null, null,null);
	}
	
	
	/**
	 * Execute a command in new thread, redirect its STDOUT and STDERR and 
	 * cancel the command if it is requested by the user via the given ExecutionContext
	 * STDOUT and STDERR are written to files
	 * STDIN is read from file
	 * @param command String array contains the program/path to the script in the first entry and all commandline parameters in the following entries
	 * @param exec ExecutionContext
	 * @param logger a NodeLogger to write messages to
	 * @param stdOutFile  File to write STDOUT to, omitted if null
	 * @param stdErrFile  File to write STDERR to, omitted if null
	 * @param sdtInFile	  File to read from STDIN, omitted if null
	 * @throws CanceledExecutionException 
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 * @throws UnsuccessfulExecutionException 
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, String [] environment, NodeLogger logger, String stdOutFile, String stdErrFile, String stdInFile) throws CanceledExecutionException, InterruptedException, ExecutionException, UnsuccessfulExecutionException{
		Executor.executeCommand(command, exec, environment, logger, stdOutFile, stdErrFile, null, null, stdInFile,null);
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
	 * @param sdtInFile	File to read from STDIN, omitted if null
	 * @throws CanceledExecutionException 
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 * @throws UnsuccessfulExecutionException 
	 * @throws Exception
	 */
	public static void executeCommand(String[] command, ExecutionContext exec, String[] environment, NodeLogger logger, String stdOutFile, String stdErrFile, StringBuffer stdOut, StringBuffer stdErr, String StdInFile,StringBuffer HTEOUT) throws CanceledExecutionException, InterruptedException, ExecutionException, UnsuccessfulExecutionException  {
		exec.checkCanceled();
		//Create Threadpool + the execution and monitoring threads 
		ExecutorService pool = Executors.newSingleThreadExecutor();
		
		ExecuteThread executorThread   = new ExecuteThread(command, logger, stdOutFile, stdErrFile, StdInFile, stdOut, stdErr, environment,HTEOUT);
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
				throw e;
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
			throw(new UnsuccessfulExecutionException("Unsuccessful Execution of " + command[0] + "\n"+IO.tail(executorThread.getSTDERR(), 5)+""));
		}
		//logger.info("ThreadPool closed successfully.");
		
		// check for exit code
		if(executorThread.getExitCode() != 0) {
			logger.error("Exit code was not 0: '"+ executorThread.getExitCode() +"'!");
			throw(new UnsuccessfulExecutionException("Exit code was not 0: '"+ executorThread.getExitCode() +"' for " + ExecuteThread.getCommand(command)));
		}
	}
	
	public static int executeCommandWithExitCode(String [] command, ExecutionContext exec, NodeLogger logger) throws CanceledExecutionException, InterruptedException, ExecutionException {
		return executeCommandWithExitCode(command, exec, null, logger, null, null, null, null,null,null);
	}
	
	/**
	 * duplicate method returning exit code
	 * @param command
	 * @param exec
	 * @param environment
	 * @param logger
	 * @param stdOutFile
	 * @param stdErrFile
	 * @param stdOut
	 * @param stdErr
	 * @param StdInFile
	 * @return
	 * @throws CanceledExecutionException
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 */
	public static int executeCommandWithExitCode(String[] command, ExecutionContext exec, String[] environment, NodeLogger logger, String stdOutFile, String stdErrFile, StringBuffer stdOut, StringBuffer stdErr, String StdInFile,StringBuffer HTEOUT) throws CanceledExecutionException, InterruptedException, ExecutionException{
		exec.checkCanceled();
		//Create Threadpool + the execution and monitoring threads 
		ExecutorService pool = Executors.newSingleThreadExecutor();
		
		ExecuteThread executorThread   = new ExecuteThread(command, logger, stdOutFile, stdErrFile, StdInFile, stdOut, stdErr, environment,HTEOUT);
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
				throw e;
			}
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				// ignore interrupted exception
			}
		}
		
		return executorThread.getExitCode();
	}
}

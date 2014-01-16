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
	    private boolean useStdOutStream;
	    private boolean useStdErrStream;
		
		/**
		 * ExecuteThread that does not! read the Error/Outstreams
		 * @param threadname
		 * @param command
		 */
		public ExecuteThread(String threadname, String command) {
		      this.name = threadname;
		      this.executingThread = new Thread(this, name);
		      System.out.println("New execute thread: " + executingThread);
		      this.finished = false;
		      this.command = command;
		      
		      this.useStdOutStream=false;
		      this.useStdErrStream=false;
		   }
		
		
		/**
		 * ExecuteThread that reads the StdOutstreams
		 * @param threadname
		 * @param command
		 * @param stdOutFile
		 */
		public ExecuteThread(String threadname, String command, String stdOutFile) {
		      this.name = threadname;
		      this.executingThread = new Thread(this, name);
		      System.out.println("New execute thread: " + executingThread);
		      this.finished = false;
		      this.command = command;
		   
		      this.stdOutFile=stdOutFile;
		      
		      this.useStdOutStream=true;
		      this.useStdErrStream=false;
		   }
		
		
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
		      
		      this.useStdOutStream=true;
		      this.useStdErrStream=true;
		   }
		
		
		   // This is the entry point for thread.
		   public void run() {
		      try {
		    	  System.out.println("Running command: "+command);
		    	  
		    	  //Start the process
		    	  p = Runtime.getRuntime().exec(command);
		    	  
		    	  //If Output is written to stdout/stderr
		    	  if(useStdOutStream){

			          StreamThread stdOutStream = new StreamThread(p.getInputStream(),stdOutFile);

			          if(useStdErrStream){//Start both streams
			        	  StreamThread stdErrStream = new StreamThread(p.getErrorStream(),stdErrFile);
			        	  stdOutStream.start();
			        	  stdErrStream.start();
				    	  //Wait for process 
						  p.waitFor();
						  //Process finished
				    	  finished = true;
				    	  //Close Stream Threads
				    	  stdOutStream.interrupt();
			        	  stdErrStream.interrupt();
			          }else{//Start Stdout only
			        	  stdOutStream.start();
						  p.waitFor();
						  //Process finished
				    	  finished = true;
				    	  //Close Stream Threads
				    	  stdOutStream.interrupt();
			          } 
		    	  }else{//No stdout/sterr threads
			    	  //Wait for process 
					  p.waitFor();
					  //Process finished
			    	  finished = true;
			    	  
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
		    
			private void getPID(Process p){
		    	try {
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

package de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutionException;

import org.apache.commons.io.FileUtils;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.port.PortType;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SettingsStorageNodeModel;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;

/**
 * Abstract class which uses Executor class to execute commands and can catch STDOUT / STDERR.
 * 
 * @author Jonas Zierer
 * @author Michel Kluge
 *
 */
public abstract class ExecutorNodeModel extends SettingsStorageNodeModel {
		
	// dummy messsages
	public static final String LOGMESSAGE_BEFORE_EXECUTION = "Not yet executed!";
	public static final String LOGMESSAGE_LOG_DISABLED     = "Not yet executed!";
	
	// StringBuffer for stdout and stderr
	private final StringBuffer STDOUT;
	private final StringBuffer STDERR;
	
	// filenames of saved stdout and stderr
	private static final String FILE_STDOUT = "stdout.log";
	private static final String FILE_STDERR = "stderr.log";

	private String stdout_custom = null;
	private String stderr_custom = null;

	// logger class
	private static final NodeLogger LOGGER = NodeLogger.getLogger(ExecutorNodeModel.class);

	//private boolean isRunning = false;
	
	/**
	 * Constructor with number of input and output ports.
	 * @param nrInDataPorts number of input ports
	 * @param nrOutDataPorts number of output ports
	 * @param catchStdout catches stdout if true
	 * @param catchStderr catches stderr if true
	 */
	protected ExecutorNodeModel(int nrInDataPorts, int nrOutDataPorts, boolean catchStdout, boolean catchStderr) {
		super(nrInDataPorts, nrOutDataPorts);
		
		if(catchStdout)
			STDOUT = new StringBuffer(LOGMESSAGE_BEFORE_EXECUTION);
		else
			STDOUT = null;
		
		if(catchStderr)
			STDERR = new StringBuffer(LOGMESSAGE_BEFORE_EXECUTION);
		else
			STDERR = null;
	}
	
	/**
	 * Constructor with number of input and output ports.
	 * @param nrInDataPorts number of input ports
	 * @param nrOutDataPorts number of output ports
	 * @param catchStdout catches stdout if true
	 * @param catchStderr catches stderr if true
	 */
	protected ExecutorNodeModel(final PortType[] inPortTypes, final PortType[] outPortTypes, boolean catchStdout, boolean catchStderr) {
		super(inPortTypes, outPortTypes);
		
		if(catchStdout)
			STDOUT = new StringBuffer(LOGMESSAGE_BEFORE_EXECUTION);
		else
			STDOUT = null;
		
		if(catchStderr)
			STDERR = new StringBuffer(LOGMESSAGE_BEFORE_EXECUTION);
		else
			STDERR = null;
	}

	/**
	 * Executes the command which is defined in command and saves stdout and stderr
	 * @param exec ExecutionContext
	 * @param command command which is joined by " " before execution
	 * @param environment environment variables
	 * @throws Exception 
	 */
	protected void executeCommand(final ExecutionContext exec, String[] command, String[] environment, File path2Stdout,  File path2Stderr) throws IOException, CanceledExecutionException, InterruptedException, ExecutionException, UnsuccessfulExecutionException  {
		//isRunning = true;
		// StringBuffers to write STDOUT and STDERR to
		if(STDOUT!=null){
			STDOUT.setLength(0);
			STDOUT.append("----------------------- START ----------------------------\n");
		}
		if(STDERR != null){
			STDERR.setLength(0);
			STDERR.append("----------------------- START ----------------------------\n");
		}
		
		// Files to write STDOUT and STDERR to
		if(path2Stdout != null) {
			if(!path2Stdout.getParentFile().exists()) {
				path2Stdout.getParentFile().mkdirs();
			}
			try {
				stdout_custom = path2Stdout.getCanonicalPath();
			} catch (IOException e) {
				LOGGER.error(e.getMessage());
				throw e;
			}
		}
		if(path2Stderr != null) {
			if(!path2Stderr.getParentFile().exists()) {
				path2Stderr.getParentFile().mkdirs();
			}
			try {
				stderr_custom = path2Stderr.getCanonicalPath();
			} catch (IOException e) {
				LOGGER.error(e.getMessage());
				throw e;
			}
		}

		Executor.executeCommand(command, exec, environment, LOGGER, stdout_custom, stderr_custom, STDOUT, STDERR);

		if(STDOUT!=null){
			STDOUT.append("------------------------ END --------------------------");
		}
		if(STDERR != null){
			STDERR.append("------------------------ END --------------------------");
		}
		//isRunning=false;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// OVERIDE KNIME NODE METHODS
	/////////////////////////////////////////////////////////////////////////////////////////////////    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
//    	if(!isRunning){ // only delete Logs if node was not interrupted during run
//			if(STDOUT!=null){
//				STDOUT.setLength(0);
//				STDOUT.append(LOGMESSAGE_BEFORE_EXECUTION);
//			}
//			
//			if(STDERR != null){
//				STDERR.setLength(0);
//				STDERR.append(LOGMESSAGE_BEFORE_EXECUTION);
//			}
//    	}
//    	isRunning = false;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	// clear the old buffers
    	this.reset();
    	// load the files if possible
    	File f_stdout = new File(internDir, FILE_STDOUT);
    	File f_stderr = new File(internDir, FILE_STDERR);
        
    	// set the buffers, if the files are there
    	if(f_stdout.exists()){
    		STDOUT.setLength(0);
    		STDOUT.append(FileUtils.readFileToString(f_stdout));
    	}
    	if(f_stderr.exists()){
    		STDERR.setLength(0);
    		STDERR.append(FileUtils.readFileToString(f_stderr));
    	}
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	// create files in the internalDir
    	File f_stdout = new File(internDir, FILE_STDOUT);
    	File f_stderr = new File(internDir, FILE_STDERR);
    	
    	// write the StringBuffers to Files
    	FileUtils.write(f_stdout, STDOUT.toString());
    	FileUtils.write(f_stderr, STDERR.toString());
    }

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// GETTERS
	/////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * 
     * @return STDOUT of the binary
     */
	public String getSTDOUT() {
		return STDOUT.toString();
	}
	
	/**
	 * 
	 * @return STDERR of the binary
	 */
	public String getSTDERR() {
		return STDERR.toString();
	}

//	protected void setSTDOUT(String stdout, boolean reset){
//		if(reset)this.STDOUT.setLength(0);
//		this.STDOUT.append(stdout);
//	}
//	protected void setSTDERR(String stderr, boolean reset){
//		if(reset)this.STDERR.setLength(0);
//		this.STDERR.append(stderr);
//	}
}

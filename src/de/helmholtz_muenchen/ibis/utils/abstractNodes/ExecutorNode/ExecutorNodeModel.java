package de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SettingsStorageNodeModel;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * Abstract class which uses Executor class to execute commands and can catch STDOUT / STDERR.
 * 
 * @author Jonas Zierer
 * @author Michel Kluge
 *
 */
public abstract class ExecutorNodeModel extends SettingsStorageNodeModel {
		
	// StringBuffer for stdout and stderr
	private final StringBuffer STDOUT;
	private final StringBuffer STDERR;
	
	// filenames of stdout and stderr
	private static final String FILE_STDOUT = "stdout.log";
	private static final String FILE_STDERR = "stderr.log";

	// logger class
	private static final NodeLogger LOGGER = NodeLogger.getLogger(ExecutorNodeModel.class);

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
			STDOUT = new StringBuffer();
		else
			STDOUT = null;
		
		if(catchStderr)
			STDERR = new StringBuffer();
		else
			STDERR = null;
	}

	/**
	 * Executes the command which is defined in command and saves stdout and stderr
	 * @param exec ExecutionContext
	 * @param command command which is joined by " " before execution
	 * @param environment environment variables
	 * @throws CanceledExecutionException
	 */
	protected void executeCommand(final ExecutionContext exec, String[] command, String[] environment) throws CanceledExecutionException {
		try {
			Executor.executeCommand(command, exec, environment, LOGGER, STDOUT, STDERR);
		} catch (Exception e) {
			LOGGER.error(e.getMessage());
			throw(new CanceledExecutionException(e.getMessage()));
		}
	}
	

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// OVERIDE KNIME NODE METHODS
	/////////////////////////////////////////////////////////////////////////////////////////////////    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
		STDOUT.delete(0, STDOUT.length());
    	STDERR.delete(0, STDERR.length());
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
    		STDOUT.append(FileUtils.readFileToString(f_stdout));
    	}
    	if(f_stderr.exists()){
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
}

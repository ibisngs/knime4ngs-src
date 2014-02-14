package de.helmholtz_muenchen.ibis.utils.abstractNodes.WrapperNode;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SettingsStorageNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode.ScriptNodeModel;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * Abstract class which can be used to implement wrapper nodes for external tools.
 * Mostly copied from ScriptNodeModel and made ScriptNodeModel extend this class in order to avoid code duplication.
 * 
 * @author Michel Kluge
 * @author Jonas Zierer
 *
 */
public abstract class WrapperNodeModel extends SettingsStorageNodeModel {
		
	// StringBuffer for stdout and stderr
	private final StringBuffer STDOUT;
	private final StringBuffer STDERR;
	
	// filenames of stdout and stderr
	private static final String FILE_STDOUT = "stdout.log";
	private static final String FILE_STDERR = "stderr.log";

	// logger class
	private static final NodeLogger LOGGER = NodeLogger.getLogger(ScriptNodeModel.class);

	/**
	 * Constructor with number of input and output ports.
	 * @param nrInDataPorts number of input ports
	 * @param nrOutDataPorts number of output ports
	 * @param catchStdout catches stdout if true
	 * @param catchStderr catches stderr if true
	 */
	protected WrapperNodeModel(int nrInDataPorts, int nrOutDataPorts, boolean catchStdout, boolean catchStderr) {
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
	
	/**
     * checks, if a binary is there and if the execution flag is set
     * @param binaryPath path to binary
     * @return
     * @throws InvalidSettingsException
     */
    protected boolean validateBinary(String binaryPath) throws InvalidSettingsException
    {
    	if(binaryPath == null || binaryPath.length() == 0)
    		throw new InvalidSettingsException("Path to binary is not set.");
    	
    	// check, if file can be executed
    	File binFile = new File(binaryPath);
    	if(!(binFile.isFile() && binFile.canExecute()))
    		throw new InvalidSettingsException("Executable flag of '" + binaryPath + "' is not set or file is not found.");
    	return true;
    }
    
    /**
     * creates a absolute path name if the path name starts with "./" or "../" which means its a relative
     * path (to the folder where the binary is located)
     * @param filenPath path to file
     * @param binaryPath File which points to the binary
     * @param isPath true, if its a path and not a filename
     * @return
     */
    protected String getAbsoluteFilename(String filenPath, File binaryPath, boolean isPath) {
    	// relative path to the binary folder
    	if(filenPath.startsWith("."))
    		filenPath = binaryPath.getParentFile().getAbsolutePath() + File.separator + filenPath;
    	
    	// check, if path ends with seperator
    	if(isPath && !filenPath.endsWith(File.separator))
    		filenPath += File.separator;
    	
    	// remove useless "/./"
    	filenPath = filenPath.replaceAll(File.separator + "\\." + File.separator, File.separator);

    	return filenPath;

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

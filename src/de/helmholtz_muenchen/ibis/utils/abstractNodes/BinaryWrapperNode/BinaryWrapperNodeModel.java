package de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode;

import java.io.File;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode.ExecutorNodeModel;


public abstract class BinaryWrapperNodeModel extends ExecutorNodeModel {

    // keys for SettingsModels
	protected static final String CFGKEY_BINARY_PATH 	= "BinaryPath";
	
	// initial default values for SettingsModels    
    private static final String DEFAULT_BINARY_PATH 	= "";
    
    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_BINARY_PATH = getSettingsModelString(CFGKEY_BINARY_PATH);
    
    
    /**
     * add the used settings
     */
    static {
        // add values for SettingsModelString
        addSettingsModelString(CFGKEY_BINARY_PATH, DEFAULT_BINARY_PATH);
    }
	
	/**
	 * Constructor with number of input and output ports.
	 * @param nrInDataPorts number of input ports
	 * @param nrOutDataPorts number of output ports
	 * @param catchStdout catches stdout if true
	 * @param catchStderr catches stderr if true
	 */
	protected BinaryWrapperNodeModel(int nrInDataPorts, int nrOutDataPorts, boolean catchStdout, boolean catchStderr) {
		super(nrInDataPorts, nrOutDataPorts, catchStdout, catchStderr);
	}
	
	@Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
		// validate binary
		validateBinary(SET_BINARY_PATH.getStringValue());
		
		return new DataTableSpec[]{null};
	}
	
	/**
	 * Returns the binary path which is currently was saved correctly
	 * @return
	 */
	protected String getBinaryPath() {
		return SET_BINARY_PATH.getStringValue();
	}
	
	/**
	 * Returns the default binary path.
	 * Can be overridden by extending if a different default value should be set.
	 * @return
	 */
	protected String getDefaultBinaryPath() {
		return DEFAULT_BINARY_PATH;
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

}

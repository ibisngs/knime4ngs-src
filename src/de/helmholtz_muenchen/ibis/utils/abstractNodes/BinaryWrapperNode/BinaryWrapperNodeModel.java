package de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModel;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.threads.ExecuteThread;

/**
 * Can be used for a node which wrapps around a binary.
 * @author Michael Kluge
 *
 */
public abstract class BinaryWrapperNodeModel extends HTExecutorNodeModel {

    // keys for SettingsModels
	protected static final String CFGKEY_BINARY_PATH 			= "BinaryPath";
	protected static final String CFGKEY_ADDITIONAL_PARAMETER 	= "AdditionalParameter";
//    protected static final String CFGKEY_PARAMETER_FILE 		= "ParameterFile";
	
	// initial default values for SettingsModels    
    protected static final String DEFAULT_BINARY_PATH 			= "";
    protected static final String DEFAULT_ADDITIONAL_PARAMETER	= "";
//    protected static final String DEFAULT_PARAMETER_FILE 			= "-"; // must be set by user but is optional
    
    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_BINARY_PATH			= new SettingsModelString(CFGKEY_BINARY_PATH, DEFAULT_BINARY_PATH);
    private final SettingsModelString SET_ADDITIONAL_PARAMETER 	= new SettingsModelString(CFGKEY_ADDITIONAL_PARAMETER, DEFAULT_ADDITIONAL_PARAMETER);
//    private final SettingsModelString SET_PARAMETER_FILE		= new SettingsModelString(CFGKEY_PARAMETER_FILE, DEFAULT_PARAMETER_FILE);
    
	// logger class
//	private static final NodeLogger LOGGER = NodeLogger.getLogger(ExecutorNodeModel.class);
	
	// regex for additional parameter
	private static final String PARAMETER_REGEX 	= "(-{1,2}\\w+)\\s*(([^-]|(?<=\\w)-)*)";
	private static final Pattern PARAMETER_PATTERN 	= Pattern.compile(PARAMETER_REGEX);
//	private static final String REGEX_WHITESPACE	= "\\s+";
	

	// stores the defined SettingsModels objects
	private final HashMap<String, SettingsModel> SETTINGS_MAP = new HashMap<String, SettingsModel>();
	private final ArrayList<SettingsModel> SETTINGS = new ArrayList<SettingsModel>();
	
	//Catch Stdout/err streams
	private boolean catchStdout;
	private boolean catchStderr;
	
	/**
	 * Constructor with number of input and output ports.
	 * @param nrInDataPorts number of input ports
	 * @param nrOutDataPorts number of output ports
	 * @param catchStdout catches stdout if true
	 * @param catchStderr catches stderr if true
	 */
	protected BinaryWrapperNodeModel(int nrInDataPorts, int nrOutDataPorts, boolean catchStdout, boolean catchStderr) {
		super(nrInDataPorts,nrOutDataPorts);
		// init is called by child classes
		
		this.catchStderr = catchStderr;
		this.catchStdout = catchStdout;
		init();
	}

	public void init() {
		super.init();
		addSetting(SET_BINARY_PATH);
		addSetting(SET_ADDITIONAL_PARAMETER);
//		addSetting(SET_PARAMETER_FILE);
	}
	
	@Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
		// validate binary
		isBinaryValid(SET_BINARY_PATH.getStringValue());
		
		return new DataTableSpec[]{null};
	}
	
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
    	
    	exec.setProgress(0.01); // tell the user that we started with the work
    	
    	// get binary path
    	String getBinaryPath = getBinaryPath();   	
    	// get Arguments
    	ArrayList<String> commands = getSetParameters(inData);
    	// add binary path as first part of the command
    	commands.add(0, getBinaryPath);
    	String[] command = commands.toArray(new String[commands.size()]);
    	System.out.println(ExecuteThread.getCommand(command));
    	// check if run was already successful 
    	File lockFile = getPathToLockFile();
    	
    	runHTExecute(command,exec,lockFile);
		
		exec.setProgress(1.00); // we are done
        return getOutputData(exec, ExecuteThread.getCommand(command), inData);
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
	 * Returns the parameters which can be set via the GUI.
	 * These will override all the other settings.
	 * @param inData Input data tables
	 * @return
	 */
	protected abstract LinkedHashMap<String, String> getGUIParameters(final BufferedDataTable[] inData);
	
	/**
	 * Returns the output data
	 * @param exec ExecutionContext
	 * @param command command which was called
	 * @return
	 */
	protected abstract BufferedDataTable[] getOutputData(final ExecutionContext exec, String command, final BufferedDataTable[] inData);
	
	/**
	 * Path to a file where STDERR should be written (or null if no log should be written)
	 * @return
	 */
	protected abstract File getPathToStderrFile();
	
	
	/**
	 * Path to a file where STDOUT should be written (or null if no log should be written)
	 * @return
	 */
	protected abstract File getPathToStdoutFile();
	
	/**
	 * Path to lock file or null if no check should be performed if node was executed sucessfully already
	 * @return
	 */
	protected abstract File getPathToLockFile();
	
	/**
	 * Returns the parameters which were set by the parameter file and the additional parameter field
	 * and by the GUI. 
	 * Override levels: GUI parameter -> additional parameter -> default parameter file
	 * @param inData Input data tables
	 * @return
	 */
	protected ArrayList<String> getSetParameters(final BufferedDataTable[] inData) {
		LinkedHashMap<String, String> pars = new LinkedHashMap<String, String>();			// merged parameter set
//		LinkedHashMap<String, String> parsFile = getParametersFromParameterFile();			// get parameter from file
		LinkedHashMap<String, String> parsAdditional = getAdditionalParameter();			// get parameter from input field
		LinkedHashMap<String, String> parsGUI = getGUIParameters(inData);					// get parameter from GUI
		
		// merge them all together
//		pars.putAll(parsFile);
		pars.putAll(parsAdditional);
		pars.putAll(parsGUI);
		
		// build the command list
		ArrayList<String> commands = new ArrayList<String>();
		for(Iterator<String> it = pars.keySet().iterator(); it.hasNext(); ) {
			// add parameter name
			String key = it.next();
			
			if(key.length() > 0) {
				// add value, if some is set
				String value = pars.get(key);
				if(value.length() != 0)
					commands.add(key + " " + value);
				else
					commands.add(key);
			}
		}
		
		// return the commands
		return commands;
	}
	
	/**
	 * splits a "complete parameter line" at the first space
	 * @param completeParameter
	 * @return
	 */
//	private Pair<String, String> splitParameter(String completeParameter) {
//		String[] split = completeParameter.split(REGEX_WHITESPACE);
//		String par = split[0];
//		String arg = ""; // default argument is empty
//		
//		if(split.length > 1)
//			arg = completeParameter.replaceFirst(split[0] + REGEX_WHITESPACE, "").replaceFirst(REGEX_WHITESPACE + "$", "");
//		return new Pair<String, String>(par, arg);
//	}
	
	/**
	 * Returns the parameter which were set in the additional parameter input field.
	 * @return
	 */
	private LinkedHashMap<String, String> getAdditionalParameter() {
		LinkedHashMap<String, String> pars = new LinkedHashMap<String, String>();
		String parameter = SET_ADDITIONAL_PARAMETER.getStringValue();
		
		// split at REGEX
		Matcher m = PARAMETER_PATTERN.matcher(parameter);
		
		while(m.find())
			pars.put(m.group(1), m.group(2));
		
		return pars;
	}
	
	/**
	 * Returns the values, which are stored in the parameter file, if one is given
	 * @return
	 */
//	private LinkedHashMap<String, String> getParametersFromParameterFile() {
//		LinkedHashMap<String, String> pars = new LinkedHashMap<String, String>();
//		String filename = SET_PARAMETER_FILE.getStringValue();
//
//		// check if some parameter file was set
//		if(!(filename.length() > 0 && !DEFAULT_PARAMETER_FILE.equals(filename)))
//			return pars;
//
//		// check, if file is there
//		File f = new File(filename);
//		if(!(f.exists() && f.isFile()))
//			return pars;
//
//		// open file
//		try {
//			String line;
//			BufferedReader r = new BufferedReader(new FileReader(f));
//			// read lines
//			while((line = r.readLine()) != null) {
//				Pair<String, String> p = splitParameter(line);
//				if(p.getFirst().length() > 0)
//					pars.put(p.getFirst(), p.getSecond());
//			}
//			// close the file
//			r.close();
//		}
//		catch(Exception e) {
//			LOGGER.error(e.getStackTrace());
//		}
//		return pars;
//	}
	
	/**
     * checks, if a binary is there and if the execution flag is set
     * @param binaryPath path to binary
     * @return
     * @throws InvalidSettingsException
     */
    protected boolean isBinaryValid(String binaryPath) throws InvalidSettingsException
    {
    	return IO.isBinaryValid(binaryPath, true);
    }
    
    /**
     * creates a absolute path name if the path name starts with "./" or "../" which means its a relative
     * path (to the folder where the binary is located)
     * @param filenPath path to file
     * @param isPath true, if its a path and not a filename
     * @return
     */
    protected String getAbsoluteFilename(String filenPath, boolean isPath) {
    	File binaryPath = new File(getBinaryPath());
     	
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
    
    
	public void addSetting(SettingsModel model) {
		if(model != null) {
			SETTINGS.add(model);
			// get protected config name
			try {
				Method m = model.getClass().getDeclaredMethod("getConfigName");
				m.setAccessible(true);
				String configName = m.invoke(model, new Object[]{}).toString();
				
				// add it to the settings map
				SETTINGS_MAP.put(configName, model);
			}
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		else {
			throw new IllegalArgumentException("SettingsModel can not be null if it should be added!");
		}
	}
    
	/**
	 * gets a saved settings or null, if not there
	 * @param name
	 */
	public SettingsModel getSetting(String name) {
		return this.SETTINGS_MAP.get(name);
	}
	
	
	
	private void runHTExecute(String[] command, ExecutionContext exec, File lockFile) throws Exception{
	
		// StringBuffer for stdout and stderr
		StringBuffer STDOUT = null;
		StringBuffer STDERR = null;
		
		if(catchStderr){
			STDERR = new StringBuffer("Not yet executed!");		
		}else if(catchStdout){
			STDOUT = new StringBuffer("Not yet executed!");		
		}
		
		
		
		//HTE Execution
		super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, null, lockFile, getPathToStdoutFile().getAbsolutePath(), getPathToStderrFile().getAbsolutePath(), STDOUT, STDERR, null);	
	}
	
	
	
	/**************************************** KNIME METHODS ****************************************/
	/***********************************************************************************************/
    
    /**
     * {@inheritDoc}
     */
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	
    	super.validateSettings(settings);
    	
    	// get all models and validate them
    	for(SettingsModel set : SETTINGS)
    		set.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	
    	super.loadValidatedSettingsFrom(settings);
    	
    	// get all models and load them
    	for(SettingsModel set : SETTINGS)
    		set.loadSettingsFrom(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	
    	super.saveSettingsTo(settings);
    	
    	// get all models and save them
    	for(SettingsModel set : SETTINGS)
    		set.saveSettingsTo(settings);
    }
	
	@Override
	protected void loadInternals(File nodeInternDir, ExecutionMonitor exec)
			throws IOException, CanceledExecutionException {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void saveInternals(File nodeInternDir, ExecutionMonitor exec)
			throws IOException, CanceledExecutionException {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void reset() {
		// TODO Auto-generated method stub
		
	}
	
}

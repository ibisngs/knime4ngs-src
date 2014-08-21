package de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelColumnName;
import org.knime.core.node.port.PortType;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode.ScriptNodeModel;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;

public abstract class RNodeModel extends ScriptNodeModel {
	public static final String ROW_ID = "ROWID"; // column name used for rowID
	
	public static final String R_SCRIPTS_PATH =  "R";
	public static final String GLOBALS_R = IO.getScriptPath() + ScriptNodeModel.SCRIPTS_SUBFOLDER + File.separatorChar + R_SCRIPTS_PATH + File.separatorChar + "utils" + File.separatorChar + "GLOBALS.R";

	private final HashMap<String,String> ARGUMENTS = new HashMap<String,String>();;
	protected final String[] INPUT_FILE_ARGUMENTS;
	protected final String[] OUTPUT_FILE_ARGUMENTS;

	private static final NodeLogger LOGGER = NodeLogger.getLogger(RNodeModel.class);
	
	protected RNodeModel(int nrInDataPorts, int nrOutDataPorts, String script, String[] input_file_arguments, String[] output_file_arguments) {
		super(nrInDataPorts, nrOutDataPorts, script);

		// check number of file names
		if(input_file_arguments.length != this.getNrInPorts()){
			LOGGER.error("HAVING " + this.getNrInPorts() + " INPUT FILES BUT "+ input_file_arguments.length + " NAMES!");
		}
		if(output_file_arguments.length != this.getNrOutPorts()){
			LOGGER.error("HAVING " + this.getNrOutPorts() + " INPUT FILES BUT "+ output_file_arguments.length + " NAMES!");
		}
		this.INPUT_FILE_ARGUMENTS = input_file_arguments;
		this.OUTPUT_FILE_ARGUMENTS= output_file_arguments;
	}
	
	protected RNodeModel(final PortType[] inPortTypes, final PortType[] outPortTypes, String script, String[] input_file_arguments, String[] output_file_arguments) {
		super(inPortTypes, outPortTypes, script);
		// TODO what if PortTypes are not tables to write???
		 
		// check number of file names
		if(input_file_arguments.length > inPortTypes.length){
			LOGGER.error("HAVING " + this.getNrInPorts() + " INPUT FILES BUT "+ input_file_arguments.length + " NAMES!");
		}
		if(output_file_arguments.length > outPortTypes.length){
			LOGGER.error("HAVING " + this.getNrOutPorts() + " INPUT FILES BUT "+ output_file_arguments.length + " NAMES!");
		}
		this.INPUT_FILE_ARGUMENTS = input_file_arguments;
		this.OUTPUT_FILE_ARGUMENTS= output_file_arguments;
	}
	

	@Override
	protected String getScriptPath(){
		return(super.getScriptPath() + R_SCRIPTS_PATH + File.separatorChar);
	}
	
	protected String getRscriptBinary(){
		if(this.getAvailableFlowVariables().containsKey("Rscript")){
			return(this.getAvailableFlowVariables().get("Rscript").getStringValue());
		}
		return("Rscript");
	}
	
	// not used in childs of this class
	@Override
	public void init() {};


	protected void prepareInputData(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{

		//////////////////////////////////////////////////////////////////////////
		// PREPARE INPUT FILES
		//////////////////////////////////////////////////////////////////////////
		exec.setProgress(0.00);
		exec.setProgress("writing input data");
		for(int i=0; i < INPUT_FILE_ARGUMENTS.length; i++){
			exec.checkCanceled();
			if(inData[i] != null){
				File tmpFile = null;
				try {
					tmpFile = File.createTempFile("knime_R_connector_" + this.SCRIPT.replaceAll(File.separatorChar+"", "_") + "_input_" , ".csv");
					this.addArgument(this.INPUT_FILE_ARGUMENTS[i], tmpFile.getCanonicalPath());
				} catch (IOException e) {
					LOGGER.error("unable to create temp file!");
					throw(e);
				}
				IO.writeAsCSV(inData[i], tmpFile, exec, LOGGER);
			}else{
				this.removeArgument(this.INPUT_FILE_ARGUMENTS[i]);
			}
			exec.checkCanceled();
		}	
	}
	
	protected String[] prepareOutputData(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{
		exec.setProgress(0.05);
		exec.setProgress("preparing output data");
		String[] outFiles = new String[OUTPUT_FILE_ARGUMENTS.length];
		for(int i=0; i < OUTPUT_FILE_ARGUMENTS.length; i ++){
			exec.checkCanceled();
			File tmpFile = null;
			try {
				tmpFile = File.createTempFile("knime_R_connector_" + this.SCRIPT.replaceAll(File.separatorChar+"", "_") + "_output_" , ".csv");
				this.addArgument(this.OUTPUT_FILE_ARGUMENTS[i], tmpFile.getCanonicalPath());
				outFiles[i] = tmpFile.getCanonicalPath();
			} catch (IOException e) {
				LOGGER.error("unable to create temp file!");
				throw(e);
			}
		}	
		return(outFiles);
	}
	
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
		//////////////////////////////////////////////////////////////////////////
		// PREPARE INPUT DATA
		//////////////////////////////////////////////////////////////////////////
		prepareInputData(inData, exec);
		
		//////////////////////////////////////////////////////////////////////////
		// PREPARE OUTPUT FILES
		//////////////////////////////////////////////////////////////////////////
		String[] outFiles = prepareOutputData(inData, exec);

		//////////////////////////////////////////////////////////////////////////
		// RUN COMMAND
		//////////////////////////////////////////////////////////////////////////
		exec.setProgress(0.10);
		exec.setProgress("executing script");
		LOGGER.info("Running Rscript with arguments: " + getArgumentsAsVector());
		try{
			super.executeScript(exec, null);
		} catch(UnsuccessfulExecutionException e){
			throw new UnsuccessfulExecutionException("R command failed\nSTDOUT:"+IO.tail(this.getSTDOUT(), 5)+"\nSTDERR:"+IO.tail(this.getSTDERR(), 5));	
		}


		//////////////////////////////////////////////////////////////////////////
		// READ DATA
		//////////////////////////////////////////////////////////////////////////
		exec.setProgress(0.90);
		exec.setProgress("reading output");
		BufferedDataTable[] output = IO.readCSV(exec, outFiles, RNodeModel.LOGGER, true, true);
		exec.setProgress(1.0);

		return(output);
	}

	public static String tail(String string){
		int i = string.lastIndexOf("\n");
		return(string.substring(i));
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// OVERIDE KNIME NODE METHODS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void reset() {
		super.reset();
		this.ARGUMENTS.clear();
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// COMMANDLINE ARGUMENTS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	public void addArgument(String arg, String value){
		this.ARGUMENTS.put(arg, value);
	}
	@SuppressWarnings("rawtypes")
	public void addArgument(String arg, Collection value){
		this.addArgument(arg, StringUtils.join(value,','));
	}
	public void addArgument(String arg, int value){
		this.addArgument(arg, String.valueOf(value));
	}
	public void addArgument(String arg, double value){
		this.addArgument(arg, String.valueOf(value));
	}
	public void addArgument(String arg, float value){
		this.addArgument(arg, String.valueOf(value));
	}
	public void addArgument(String arg, SettingsModelColumnName value){
		if(value.getColumnName() == null){
			if(value.useRowID()){
				this.addArgument(arg, RNodeModel.ROW_ID);	
			}
		}else{
			this.addArgument(arg, value.getColumnName());
		}
	}
	public void removeArgument(String arg){
		this.ARGUMENTS.remove(arg);
	}

	public void addFlag(String arg){
		this.ARGUMENTS.put(arg, null);
	}
	public void removeFlag(String arg){
		this.removeArgument(arg);
	}
	public void setFlag(String arg, boolean b){
		if(b){
			this.addFlag(arg);
		}else{
			this.removeFlag(arg);	
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// GETTERS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	public String getArguments(){
		StringBuilder args = new StringBuilder(" --globals \"" + GLOBALS_R + "\"");
		for (String key : this.ARGUMENTS.keySet()) {
			// key
			args.append(" " + key);

			// value
			if(this.ARGUMENTS.get(key) != null){
				args.append(" \"" + this.ARGUMENTS.get(key) + "\"");
			}
		}
		return(args.toString());
	}

	@Override
	protected String[] getCommand(){
		ArrayList<String> command = new ArrayList<String>();
		command.add(this.getRscriptBinary());
		command.add(this.SCRIPT);
		command.add("--globals");
		command.add(GLOBALS_R);
		for (String key : this.ARGUMENTS.keySet()) {
			command.add(key);
			// value
			if(this.ARGUMENTS.get(key) != null){
				command.add(this.ARGUMENTS.get(key));
			}
		}

		return(command.toArray(new String[command.size()]));
	}


	public String getArgumentsAsVector(){
		StringBuilder args = new StringBuilder("args.in <- c( \"--globals\", \"" + GLOBALS_R + "\"");
		for (String key : this.ARGUMENTS.keySet()) {
			// key
			args.append(", \"" + key + "\"");

			// value
			if(this.ARGUMENTS.get(key) != null){
				args.append(", \"" + this.ARGUMENTS.get(key) + "\"");
			}
		}
		args.append(")");
		return(args.toString());
	}



}

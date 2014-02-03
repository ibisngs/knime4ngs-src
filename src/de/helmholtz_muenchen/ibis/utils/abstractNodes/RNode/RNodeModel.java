package de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode.ScriptNodeModel;

public abstract class RNodeModel extends ScriptNodeModel {

	public static final String R_SCRIPTS_PATH =  "scripts" + File.separatorChar + "R" + File.separatorChar;
	public static final String GLOBALS_R = IO.getScriptPath() + R_SCRIPTS_PATH + "utils" + File.separatorChar + "GLOBALS.R";;

	private final HashMap<String,String> ARGUMENTS;
	protected final String[] INPUT_FILE_ARGUMENTS;
	protected final String[] OUTPUT_FILE_ARGUMENTS;

	private static final NodeLogger LOGGER = NodeLogger.getLogger(RNodeModel.class);

	protected RNodeModel(int nrInDataPorts, int nrOutDataPorts, String script, String[] input_file_arguments, String[] output_file_arguments) {
		super(nrInDataPorts, nrOutDataPorts, script);

		// init arguments
		this.ARGUMENTS = new HashMap<String,String>();

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

	@Override
	protected String getScriptPath(){
		return(IO.getScriptPath() + R_SCRIPTS_PATH );
	}
	
	protected String getRscriptBinary(){
		if(!this.getAvailableFlowVariables().containsKey("Rscript")){
			this.pushFlowVariableString("Rscript", "Rscript");
		}
		return(this.getAvailableFlowVariables().get("Rscript").getStringValue());
	}


	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
		
		//////////////////////////////////////////////////////////////////////////
		// PREPARE INPUT FILES
		//////////////////////////////////////////////////////////////////////////
		exec.setProgress(0.00);
		exec.setProgress("writing input data");
		for(int i=0; i < this.getNrInPorts(); i++){
			File tmpFile = null;
			try {
				tmpFile = File.createTempFile("knime_R_connector_" + this.SCRIPT.replaceAll(File.separatorChar+"", "_") + "_input_" , ".csv");
				this.addArgument(this.INPUT_FILE_ARGUMENTS[i], tmpFile.getCanonicalPath());
			} catch (IOException e) {
				LOGGER.error("unable to create temp file!");
				throw(new CanceledExecutionException("unable to create temp file!" + e.getMessage()));
			}
			IO.writeAsCSV(inData[i], tmpFile, exec, LOGGER);
		}	

		//////////////////////////////////////////////////////////////////////////
		// PREPARE OUTPUT FILES
		//////////////////////////////////////////////////////////////////////////
		exec.setProgress(0.05);
		exec.setProgress("preparing output data");
		String[] outFiles = new String[this.getNrOutPorts()];
		for(int i=0; i < this.getNrOutPorts(); i ++){
			File tmpFile = null;
			try {
				tmpFile = File.createTempFile("knime_R_connector_" + this.SCRIPT.replaceAll(File.separatorChar+"", "_") + "_output_" , ".csv");
				this.addArgument(this.OUTPUT_FILE_ARGUMENTS[i], tmpFile.getCanonicalPath());
				outFiles[i] = tmpFile.getCanonicalPath();
			} catch (IOException e) {
				LOGGER.error("unable to create temp file!");
				throw(new CanceledExecutionException("unable to create temp file!" + e.getMessage()));
			}
		}	

		//////////////////////////////////////////////////////////////////////////
		// RUN COMMAND
		//////////////////////////////////////////////////////////////////////////
		exec.setProgress(0.10);
		exec.setProgress("executing script");
		LOGGER.info("Running Rscript with arguments: " + getArgumentsAsVector());
		super.executeScript(exec, null);

		//////////////////////////////////////////////////////////////////////////
		// READ DATA
		//////////////////////////////////////////////////////////////////////////
		exec.setProgress(0.90);
		exec.setProgress("reading output");
		BufferedDataTable[] output = IO.readCSV(exec, outFiles, RNodeModel.LOGGER, true, true);
		exec.setProgress(1.0);

		return(output);
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

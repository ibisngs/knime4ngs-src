package de.helmholtz_muenchen.ibis.ngs.star;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.ngs.fastaSelector.FastaSelectorNodeModel;
import de.helmholtz_muenchen.ibis.ngs.runaligner.RunAlignerNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.WrapperNode.WrapperNodeModel;


/**
 * This is the model implementation for the wrapper of STAR.
 * STAR aligns RNA-seq reads to a reference genome using uncompressed suffix arrays. 
 * For details, please see paper: 
 * Dobin et al, Bioinformatics 2012; doi: 10.1093/bioinformatics/bts635
 *
 * @author Michael Kluge
 */
public class StarNodeModel extends WrapperNodeModel {
     
	// name of the output variables
	public static final String OUTPUT_NAME_RUN_MODE = "RunMode";
	public static final String OUTPUT_NAME_OUTPUT_FOLDER = "OutputFolder";
	public static final String OUTPUT_NAME_PARAMETER_CONFIG_FILE = "ParameterConfigFile";

    // keys for SettingsModels
    protected static final String CFGKEY_BINARY_PATH 	= "BinaryPath";
    protected static final String CFGKEY_PARAMETER_FILE = "ParameterFile";
    protected static final String CFGKEY_RUN_MODE 		= "RunMode";
    protected static final String CFGKEY_OUTPUT_FOLDER 	= "OutputFolder";
    protected static final String CFGKEY_GENOME_FOLDER 	= "GenomeFolder";

    // initial default values for SettingsModels    
    protected static final String DEFAULT_RUN_MODE 		= "alignReads";			// align read mode is default
    protected static final String ALTERNATIVE_RUN_MODE 	= "genomeGenerate";		// alternative run mode
    private static final String DEFAULT_BINARY_PATH 	= "";
    private static final String DEFAULT_PARAMETER_FILE 	= "-";					// must be set by user but is optional
    private static final String DEFAULT_OUTPUT_FOLDER 	= "./output/";			// creates a folder "output" relative to the STAR binary
    private static final String DEFAULT_GENOME_FOLDER 	= "./GenomeDir/";		// searches for the genome index in "GenomeDir" relative to the STAR binary

    // name of parameters which are defined in the STAR binary
    private final static String NAME_OF_PARAMETER_FILE_PARAM 	= "--parametersFiles";		// parameter in STAR which allows the user to define a parameter File
    private final static String NAME_OF_OUTPUT_PREFIX_PARAM 	= "--outFileNamePrefix";	// parameter in STAR which allows the user to set the output folder (relative or absolute)
    private final static String NAME_OF_OUTPUT_GENOMEDIR_PARAM 	= "--genomeDir";			// output parameter in case of generateGenome run or input in other case
    private final static String NAME_OF_RUN_MODE_PARAM 			= "--runMode";				// parameter in STAR which sets the runMode
    private final static String NAME_OF_FASTAQ_PARAM 			= "--readFilesIn";			// string(s): paths to files that contain input read1 (and, if needed,  read2) 
	
    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_BINARY_PATH;
    private final SettingsModelString SET_PARAMETER_FILE;
    private final SettingsModelString SET_RUN_MODE;
    private final SettingsModelString SET_OUTPUT_FOLDER;
    private final SettingsModelString SET_GENOME_FOLDER;
    
    // the logger instance
    @SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(StarNodeModel.class);
       
    /**
     * Constructor for the node model.
     */
    protected StarNodeModel() {
        super(1, 1, true, true);
        
        // add values for SettingsModelString
        addSettingsModelString(CFGKEY_BINARY_PATH, DEFAULT_BINARY_PATH);
        addSettingsModelString(CFGKEY_PARAMETER_FILE, DEFAULT_PARAMETER_FILE);
        addSettingsModelString(CFGKEY_RUN_MODE, DEFAULT_RUN_MODE);
        addSettingsModelString(CFGKEY_OUTPUT_FOLDER, DEFAULT_OUTPUT_FOLDER);
        addSettingsModelString(CFGKEY_GENOME_FOLDER, DEFAULT_GENOME_FOLDER);
        
        // assign the SettingsModel
        SET_BINARY_PATH 	= getSettingsModelString(CFGKEY_BINARY_PATH);
        SET_PARAMETER_FILE 	= getSettingsModelString(CFGKEY_PARAMETER_FILE);
        SET_RUN_MODE 		= getSettingsModelString(CFGKEY_RUN_MODE);
        SET_OUTPUT_FOLDER 	= getSettingsModelString(CFGKEY_OUTPUT_FOLDER);
        SET_GENOME_FOLDER 	= getSettingsModelString(CFGKEY_GENOME_FOLDER);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	// check input port
    	String[] cn=inSpecs[0].getColumnNames();
    	if(isAlignRunMode()) {
    		if(!(cn[0].equals(RunAlignerNodeModel.OUTPUT_NAME_READ_FILE1) && cn[1].equals(RunAlignerNodeModel.OUTPUT_NAME_READ_FILE2)))
    			throw new InvalidSettingsException("Incompatible input: In 'alignReads' mode the node expects the output of a 'RunAligner' node.");
    	}
    	else {
    		if(!cn[0].equals(FastaSelectorNodeModel.OUTPUT_NAME_FASTA_FILES))
    			throw new InvalidSettingsException("Incompatible input: In 'genomeGenerate' mode the node expects the output of a 'FastaSelector' node.");
    	}
    	
    	validateBinary(SET_BINARY_PATH.getStringValue());
    	
		return new DataTableSpec[]{null};
    }
    

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
    	
    	exec.setProgress(0.01); // start with the work
    	// some simple arguments
    	String binaryPath = SET_BINARY_PATH.getStringValue();
    	File binaryFile = new File(binaryPath);
    	ArrayList<String> runmodeArgument = new ArrayList<String>();
    	runmodeArgument.add(NAME_OF_RUN_MODE_PARAM);
    	runmodeArgument.add(SET_RUN_MODE.getStringValue());
    	
    	/********************* OUTPUT ****************************/
    	ArrayList<String> outputArgument = new ArrayList<String>();
    	// get output folder info and create the folder if needed
    	String outputFolder = getAbsoluteFilename(SET_OUTPUT_FOLDER.getStringValue(), binaryFile, true);
    	File outDir = new File(outputFolder);
    	if(!outDir.isDirectory())
    		outDir.mkdirs();
    	
    	
    	// check, which kind of output parameter must be set.
    	String outputFolderParameter = NAME_OF_OUTPUT_PREFIX_PARAM;
    	if(!isAlignRunMode())
    		outputFolderParameter = NAME_OF_OUTPUT_GENOMEDIR_PARAM;
    	
    	outputArgument.add(outputFolderParameter);
    	outputArgument.add(outputFolder);
    	/*********************************************************/
    	
    	/********************** INPUT ****************************/
    	ArrayList<String> inputArgument = new ArrayList<String>();
    	inputArgument.add(NAME_OF_FASTAQ_PARAM);   	
    	String parameterFile = "";
    	
    	// get input parameter from runaligner
    	if(isAlignRunMode()) {
	    	String path2readFile1 = inData[0].iterator().next().getCell(0).toString();
	    	String path2readFile2 = inData[0].iterator().next().getCell(1).toString();
	    	// add first input file
	    	inputArgument.add(getAbsoluteFilename(path2readFile1, binaryFile, false));
	    	
	    	// add second input file, if paired mode was selected and both files are different
	    	if(!path2readFile1.equals(path2readFile2) && path2readFile2.length() > 0)
	    		inputArgument.add(getAbsoluteFilename(path2readFile2, binaryFile, false));
	    	
	    	// add genome Folder to parameters
	    	inputArgument.add(NAME_OF_OUTPUT_GENOMEDIR_PARAM);
	    	inputArgument.add(getAbsoluteFilename(SET_GENOME_FOLDER.getStringValue(), binaryFile, true));
	    	
	    	// optional parameter file
	    	if(!DEFAULT_PARAMETER_FILE.equals(SET_PARAMETER_FILE.getStringValue())) {	  
	    		parameterFile = getAbsoluteFilename(SET_PARAMETER_FILE.getStringValue(), binaryFile, false);
	    	    inputArgument.add(NAME_OF_PARAMETER_FILE_PARAM);
		    	inputArgument.add(parameterFile);
	    	}
	    	    	    	
    	}
    	// get input parameter from FastaSelector (which are already absolute)
    	else {
    		for(Iterator<DataRow> it = inData[0].iterator(); it.hasNext(); )
    			inputArgument.add(it.next().getCell(0).toString());
    	}
    	/*********************************************************/
    	
    	/******************* RUN COMMAND *************************/
    	// collect all commands
    	ArrayList<String> allCommands = new ArrayList<String>();
    	allCommands.add(binaryPath);
    	allCommands.addAll(runmodeArgument);
    	allCommands.addAll(inputArgument);
    	allCommands.addAll(outputArgument);
    	
		// build command array
		String[] command = allCommands.toArray(new String[allCommands.size()]);
		
		// execute the command
		exec.setProgress(0.05);
		executeCommand(exec, command, null);
		exec.setProgress(1.00); // we are done
		/*********************************************************/

		/******************* PREPARE OUTPUT **********************/
		DataColumnSpecCreator col1 = new DataColumnSpecCreator(OUTPUT_NAME_RUN_MODE, StringCell.TYPE);
		DataColumnSpecCreator col2 = new DataColumnSpecCreator(OUTPUT_NAME_OUTPUT_FOLDER, StringCell.TYPE);
		DataColumnSpecCreator col3 = new DataColumnSpecCreator(OUTPUT_NAME_PARAMETER_CONFIG_FILE, StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec(), col2.createSpec(), col3.createSpec()};
    	DataTableSpec table = new DataTableSpec(cols);
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	StringCell cl1 = new StringCell(SET_RUN_MODE.getStringValue());
    	StringCell cl2 = new StringCell(outputFolder);
    	StringCell cl3 = new StringCell(parameterFile);
    	DataCell[] c = new DataCell[]{cl1, cl2, cl3};
    	DefaultRow r = new DefaultRow("Row0", c);
    	cont.addRowToTable(r);
    	cont.close();
    	BufferedDataTable out = cont.getTable();
		
        return new BufferedDataTable[]{out};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	
        SET_BINARY_PATH.saveSettingsTo(settings);
        SET_PARAMETER_FILE.saveSettingsTo(settings);
        SET_RUN_MODE.saveSettingsTo(settings);
        SET_OUTPUT_FOLDER.saveSettingsTo(settings);
        SET_GENOME_FOLDER.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	
        SET_BINARY_PATH.loadSettingsFrom(settings);
        SET_PARAMETER_FILE.loadSettingsFrom(settings);
        SET_RUN_MODE.loadSettingsFrom(settings);
        SET_OUTPUT_FOLDER.loadSettingsFrom(settings);
        SET_GENOME_FOLDER.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
     	
        SET_BINARY_PATH.validateSettings(settings);
        SET_PARAMETER_FILE.validateSettings(settings);
        SET_RUN_MODE.validateSettings(settings);
        SET_OUTPUT_FOLDER.validateSettings(settings);
        SET_GENOME_FOLDER.validateSettings(settings);
        
        // get the value even if it is not saved and validate it.
        validateBinary(((SettingsModelString) SET_BINARY_PATH.createCloneWithValidatedValue(settings)).getStringValue());
    }
    
    /**
     * true, if align run mode is configured
     * @return
     */
    private boolean isAlignRunMode() {
    	return DEFAULT_RUN_MODE.equals(SET_RUN_MODE.getStringValue());
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected boolean validateBinary(String binaryPath) throws InvalidSettingsException
    {
    	boolean ret = super.validateBinary(binaryPath);
    	
    	// TODO: check, if this is really a STAR binary, but HOW ?!
    	return ret;
    }
}


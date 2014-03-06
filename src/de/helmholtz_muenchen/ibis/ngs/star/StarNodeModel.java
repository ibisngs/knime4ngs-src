package de.helmholtz_muenchen.ibis.ngs.star;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;

import org.apache.commons.lang3.StringUtils;
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

import de.helmholtz_muenchen.ibis.ngs.fastaSelector.FastaSelectorNodeModel;
import de.helmholtz_muenchen.ibis.ngs.rawreadmanipulator.RawReadManipulatorNodeModel;
import de.helmholtz_muenchen.ibis.ngs.runaligner.RunAlignerNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeModel;

/**
 * This is the model implementation for the wrapper of STAR.
 * STAR aligns RNA-seq reads to a reference genome using uncompressed suffix arrays. 
 * For details, please see paper: 
 * Dobin et al, Bioinformatics 2012; doi: 10.1093/bioinformatics/bts635
 *
 * @author Michael Kluge
 */
public class StarNodeModel extends BinaryWrapperNodeModel {
     
	// name of the output variables
	public static final String OUT_COL1 = "RunMode";
	public static final String OUT_COL2 = "OutputFolder";
	public static final String OUT_COL3 = "CallCommand";

    // keys for SettingsModels
    protected static final String CFGKEY_RUN_MODE 		= "RunMode";
    protected static final String CFGKEY_OUTPUT_FOLDER 	= "OutputFolder";
    protected static final String CFGKEY_GENOME_FOLDER 	= "GenomeFolder";

    // initial default values for SettingsModels    
    protected static final String DEFAULT_RUN_MODE 		= "alignReads";			// align read mode is default
    protected static final String ALTERNATIVE_RUN_MODE 	= "genomeGenerate";		// alternative run mode
    private static final String DEFAULT_OUTPUT_FOLDER 	= "./output/";			// creates a folder "output" relative to the STAR binary
    private static final String DEFAULT_GENOME_FOLDER 	= "./GenomeDir/";		// searches for the genome index in "GenomeDir" relative to the STAR binary

    // name of parameters which are defined in the STAR binary
    private final static String NAME_OF_OUTPUT_PREFIX_PARAM 	= "--outFileNamePrefix";	// parameter in STAR which allows the user to set the output folder (relative or absolute)
    private final static String NAME_OF_OUTPUT_GENOMEDIR_PARAM 	= "--genomeDir";			// output parameter in case of generateGenome run or input in other case
    private final static String NAME_OF_RUN_MODE_PARAM 			= "--runMode";				// parameter in STAR which sets the runMode
    private final static String NAME_OF_FASTAQ_PARAM 			= "--readFilesIn";			// string(s): paths to files that contain input read1 (and, if needed,  read2)
    private final static String NAME_OF_FASTA_FILES_PARAM		= "--genomeFastaFiles";		// fasta files with genomic sequence for genome files generation
    
    private final static String NAME_OF_GENOME_PARAMETER_FILE	= "genomeParameters.txt"; 	// name of settings file from a indexed genome 
	
    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_RUN_MODE			= getSettingsModelString(CFGKEY_RUN_MODE);
    private final SettingsModelString SET_OUTPUT_FOLDER		= getSettingsModelString(CFGKEY_OUTPUT_FOLDER);
    private final SettingsModelString SET_GENOME_FOLDER		= getSettingsModelString(CFGKEY_GENOME_FOLDER);
    
    // the logger instance
    @SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(StarNodeModel.class);
       
    /**
     * add the used settings
     */
    static {
        // add values for SettingsModelString
        addSettingsModelString(CFGKEY_RUN_MODE, DEFAULT_RUN_MODE);
        addSettingsModelString(CFGKEY_OUTPUT_FOLDER, DEFAULT_OUTPUT_FOLDER);
        addSettingsModelString(CFGKEY_GENOME_FOLDER, DEFAULT_GENOME_FOLDER);
    }
    
    /**
     * Constructor for the node model.
     */
    protected StarNodeModel() {
        super(1, 1, true, true);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	super.configure(inSpecs);
    	
    	// check input port
    	String[] cn=inSpecs[0].getColumnNames();
    	if(isAlignRunMode()) {
    		if(!((cn[0].equals(RunAlignerNodeModel.OUT_COL1) && cn[1].equals(RunAlignerNodeModel.OUT_COL2) ||
    			(cn[0].equals(RawReadManipulatorNodeModel.OUT_COL1) && cn[1].equals(RawReadManipulatorNodeModel.OUT_COL2)))))
    			throw new InvalidSettingsException("Incompatible input: In 'alignReads' mode the node expects the output of a 'RunAligner' or 'RawReadManipulator' node.");
    		
            // validate genome dir
            validateGenomeIndex(SET_GENOME_FOLDER.getStringValue());
    	}
    	else {
    		if(!cn[0].equals(FastaSelectorNodeModel.OUTPUT_NAME_FASTA_FILES))
    			throw new InvalidSettingsException("Incompatible input: In 'genomeGenerate' mode the node expects the output of a 'FastaSelector' node.");
    	}
    	
		return new DataTableSpec[]{null};
    }
    
    

	@Override
	protected LinkedHashMap<String, String> getGUIParameters(final BufferedDataTable[] inData) {
		LinkedHashMap<String, String> pars = new LinkedHashMap<String, String>();

		/********************* RUN MODE ***************************/
		pars.put(NAME_OF_RUN_MODE_PARAM, SET_RUN_MODE.getStringValue());
		
		
		/********************* OUTPUT ****************************/
    	// check, which kind of output parameter must be set.
    	String outputFolderParameter = (isAlignRunMode() ? NAME_OF_OUTPUT_PREFIX_PARAM : NAME_OF_OUTPUT_GENOMEDIR_PARAM);
    	
		String outputFolderArgument = getAbsoluteFilename(SET_OUTPUT_FOLDER.getStringValue(), true);
    	File outDir = new File(outputFolderArgument);
    	// create folder, if not already there
    	if(!outDir.isDirectory())
    		outDir.mkdirs();
    	
    	pars.put(outputFolderParameter, outputFolderArgument);

    	
    	/********************** INPUT ****************************/
    	String inputParameter = (isAlignRunMode() ? NAME_OF_FASTAQ_PARAM : NAME_OF_FASTA_FILES_PARAM);
    	ArrayList<String> inputArgument = new ArrayList<String>();
    	
    	// get input parameter from run aligner
    	if(isAlignRunMode()) {	
	    	String path2readFile1 = inData[0].iterator().next().getCell(0).toString();
	    	String path2readFile2 = inData[0].iterator().next().getCell(1).toString();
	    	// add first input file
	    	inputArgument.add(getAbsoluteFilename(path2readFile1, false));
	    	
	    	// add second input file, if paired mode was selected and both files are different
	    	if(!path2readFile1.equals(path2readFile2) && path2readFile2.length() > 0)
	    		inputArgument.add(getAbsoluteFilename(path2readFile2, false));
	    	
	    	// add genome folder to parameters
	    	pars.put(NAME_OF_OUTPUT_GENOMEDIR_PARAM, SET_GENOME_FOLDER.getStringValue());   	    	
    	}
    	// get input parameter from FastaSelector (which are already absolute)
    	else {
    		for(Iterator<DataRow> it = inData[0].iterator(); it.hasNext(); )
    			inputArgument.add(it.next().getCell(0).toString());
    	}
    	// add the input parameter
    	pars.put(inputParameter, StringUtils.join(inputArgument, " "));
    	// return the GUI parameter
		return pars;
	}
	
	
	@Override
	protected BufferedDataTable[] getOutputData(final ExecutionContext exec, String command, final BufferedDataTable[] inData) {
		BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, StringCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, StringCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL3, StringCell.TYPE).createSpec()}));
    	
    	DataCell[] c = new DataCell[]{
    			new StringCell(SET_RUN_MODE.getStringValue()),
    			new StringCell(getAbsoluteFilename(SET_OUTPUT_FOLDER.getStringValue(), true)),
    			new StringCell(command)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
		
        return new BufferedDataTable[]{cont.getTable()};
	}
    
    /**
     * true, if align run mode is configured
     * @return
     */
    private boolean isAlignRunMode() {
    	return DEFAULT_RUN_MODE.equals(SET_RUN_MODE.getStringValue());
    }
    
    
    /**
     * Checks, if the path contains a genome index for STAR
     * @param genomePath
     * @return
     * @throws InvalidSettingsException
     */
    protected boolean validateGenomeIndex(String genomePath) throws InvalidSettingsException {
    	// check for relative path if binary file is valid
    	if(isBinaryValid(getBinaryPath()))
    		genomePath = getAbsoluteFilename(genomePath, true);
    	
    	// check if genomePath exists
    	File f = new File(genomePath);
    	if(!(f.isDirectory() && f.exists()))
    		throw new InvalidSettingsException("Genome path '" + genomePath + "' does not exist.");
    	
    	// check if parameter file is there
    	f = new File(genomePath + File.separator + NAME_OF_GENOME_PARAMETER_FILE);
    	if(!(f.isFile() && f.canRead()))
    		throw new InvalidSettingsException("Folder was found but it seems that it does not contain a valid genome index.");    	
    	// all checks where ok
    	return true;
    }
}


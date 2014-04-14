package de.helmholtz_muenchen.ibis.ngs.featureCounts;

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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeModel;

/**
 * This is the model implementation for the wrapper of featureCounts.
 * featureCounts: an efficient general purpose program for assigning sequence reads to genomic features
 * For details, please see paper: 
 * Liao et al, Bioinformatics 2013; doi: 10.1093/bioinformatics/btt656
 *
 * @author Michael Kluge
 */
public class FeatureCountsNodeModel extends BinaryWrapperNodeModel {
     
	// name of the output variables
	public static final String OUT_COL1 = "OutputFile";
	public static final String OUT_COL2 = "CallCommand";

    // keys for SettingsModels
    protected static final String CFGKEY_OUTPUT_FILE 				= "OutputFolder";
    protected static final String CFGKEY_ANNOTATION_FILE			= "AnnotationFile";
    protected static final String CFGKEY_ANNOTATION_TYPE			= "AnnotationType";
    protected static final String CFGKEY_ANNOTATION_FEATURE			= "AnnotationFeature";
    protected static final String CFGKEY_THREAD_NUMBER				= "ThreadNumber";
    protected static final String CFGKEY_COUNT_MULTIMAPPED			= "MultimappedFlag";
    protected static final String CFGKEY_COUNT_OVERLAPPING_MULTI	= "MultiOverlapping";
    protected static final String CFGKEY_COUNT_FRAGMENTS			= "Fragments";
    protected static final String CFGKEY_COUNT_CHIMERIC_FRAGMENTS	= "ChimericFragments";
    
    // initial default values for SettingsModels
    protected static final String DEFAULT_ANNOTATION_TYPE			= "GTF";			// input of GTF file as annotation is default value
    protected static final String ALTERNATIVE_ANNOTATION_TYPE 		= "SAF";			// alternative annotation type
    private static final String DEFAULT_OUTPUT_FOLDER 				= "./output/";		// creates a folder "output" relative to the STAR binary
    private static final String DEFAULT_ANNOTATION_FILE				= "";
    private static final String DEFAULT_ANNOTATION_FEATURE			= "exon";			// default feature which is used for counting
    private static final int DEFAULT_THREAD_NUMBER					= 1;				// default threads to use
    private static final boolean DEAFULT_COUNT_MULTIMAPPED			= false;			// do not count multimapped reads by default
    private static final boolean DEAFULT_COUNT_MULTI_OVERLAPING		= false;			// do not count multi overlapping reads by default
    private static final boolean DEAFULT_COUNT_FRAGMENTS			= false;			// only for paired reads
    private static final boolean DEAFULT_COUNT_CHIMERIC_FRAGMENTS 	= false;			// do not count chimeric fragments
    
    // name of parameters which are defined in the featureCounts binary
    private final static String NAME_OF_OUTPUT_FILE 		= "-o";		// output file
    private final static String NAME_OF_ANNOTATION_FILE 	= "-a";		// input annotation file
    private final static String NAME_OF_THREAD_NUMBER		= "-T";		// sets the number of threads
    private final static String NAME_OF_FEATURE_TYPE		= "-t";		// sets the name of the feature which is used for counting
    private final static String NAME_OF_COUNT_MULTIMAPPED	= "-M";		// sets if multimapped reads should be counted or not
    private final static String NAME_OF_ANNOTATION_TYPE		= "-F";		// annotation type parameter
    private final static String NAME_OF_ASSIGN_MULTI		= "-O";		// enables counting for more than one fragment
    private final static String NAME_OF_COUNT_FRAGMENTS		= "-p";		// enables count of fragments
    private final static String NAME_OF_COUNT_CHIMERIC		= "-C";		// disables count of chimeric fragments
    				
    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_FEATURE_TYPE				= getSettingsModelString(CFGKEY_ANNOTATION_FEATURE);
    private final SettingsModelString SET_OUTPUT_FILE 				= getSettingsModelString(CFGKEY_OUTPUT_FILE);
    private final SettingsModelString SET_ANNOTATION_FILE			= getSettingsModelString(CFGKEY_ANNOTATION_FILE);
    private final SettingsModelString SET_ANNOTATION_TYPE			= getSettingsModelString(CFGKEY_ANNOTATION_TYPE);
    private final SettingsModelInteger SET_THREAD_NUMBER			= getSettingsModelInteger(CFGKEY_THREAD_NUMBER);
    private final SettingsModelBoolean SET_COUNT_MULTIMAPPED		= getSettingsModelBoolean(CFGKEY_COUNT_MULTIMAPPED);
    private final SettingsModelBoolean SET_COUNT_OVERLAPPING_MULTI	= getSettingsModelBoolean(CFGKEY_COUNT_OVERLAPPING_MULTI);
    private final SettingsModelBoolean SET_COUNT_FRAGMENTS			= getSettingsModelBoolean(CFGKEY_COUNT_FRAGMENTS);
    private final SettingsModelBoolean SET_CHIMERIC_FRAGMENTS		= getSettingsModelBoolean(CFGKEY_COUNT_CHIMERIC_FRAGMENTS);

    protected final static int MIN_THREADS = 1;
    protected final static int MAX_THREADS = 16;
    
    // the logger instance
    @SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FeatureCountsNodeModel.class);
       
    /**
     * add the used settings
     */
    static {
        // add values for SettingsModelString
    	addSettingsModelString(CFGKEY_OUTPUT_FILE, DEFAULT_OUTPUT_FOLDER);
        addSettingsModelString(CFGKEY_ANNOTATION_TYPE, DEFAULT_ANNOTATION_TYPE);
        addSettingsModelString(CFGKEY_ANNOTATION_FILE, DEFAULT_ANNOTATION_FILE);
        addSettingsModelInteger(CFGKEY_THREAD_NUMBER, DEFAULT_THREAD_NUMBER);
        addSettingsModelBoolean(CFGKEY_COUNT_MULTIMAPPED, DEAFULT_COUNT_MULTIMAPPED);
        addSettingsModelBoolean(CFGKEY_COUNT_OVERLAPPING_MULTI, DEAFULT_COUNT_MULTI_OVERLAPING);
        addSettingsModelBoolean(CFGKEY_COUNT_FRAGMENTS, DEAFULT_COUNT_FRAGMENTS);
        addSettingsModelBoolean(CFGKEY_COUNT_CHIMERIC_FRAGMENTS, DEAFULT_COUNT_CHIMERIC_FRAGMENTS);
        addSettingsModelString(CFGKEY_ANNOTATION_FEATURE, DEFAULT_ANNOTATION_FEATURE);
    }
    
    /**
     * Constructor for the node model.
     */
    protected FeatureCountsNodeModel() {
        super(1, 1, true, true);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	super.configure(inSpecs);
    	
        validateAnnotationFile(SET_ANNOTATION_FILE.getStringValue());

		return new DataTableSpec[]{null};
    }
    
    

	@Override
	protected LinkedHashMap<String, String> getGUIParameters(final BufferedDataTable[] inData) {
		LinkedHashMap<String, String> pars = new LinkedHashMap<String, String>();

		/********************* SIMPLE PARAMETER ***************************/
		if(SET_COUNT_MULTIMAPPED.getBooleanValue())
			pars.put(NAME_OF_COUNT_MULTIMAPPED, "");
		if(SET_COUNT_OVERLAPPING_MULTI.getBooleanValue())
			pars.put(NAME_OF_ASSIGN_MULTI, "");
		if(SET_COUNT_FRAGMENTS.getBooleanValue())
			pars.put(NAME_OF_COUNT_FRAGMENTS, "");
		if(!SET_CHIMERIC_FRAGMENTS.getBooleanValue())
			pars.put(NAME_OF_COUNT_CHIMERIC, "");
		if(SET_THREAD_NUMBER.getIntValue() > 1)
			pars.put(NAME_OF_THREAD_NUMBER, Integer.toString(SET_THREAD_NUMBER.getIntValue()));
		
		pars.put(NAME_OF_ANNOTATION_TYPE, SET_ANNOTATION_TYPE.getStringValue());
		pars.put(NAME_OF_FEATURE_TYPE, SET_FEATURE_TYPE.getStringValue());
		pars.put(NAME_OF_ANNOTATION_FILE, SET_ANNOTATION_FILE.getStringValue());
		
		/********************* OUTPUT ****************************/
		String outputFolderArgument = getAbsoluteFilename(SET_OUTPUT_FILE.getStringValue(), false);
    	File outDir = new File(outputFolderArgument).getParentFile();
    	// create folder, if not already there
    	if(!outDir.isDirectory())
    		outDir.mkdir();
    	
    	pars.put(NAME_OF_OUTPUT_FILE, outputFolderArgument);

    	
    	/********************** INPUT BAM/SAM ****************************/
    	ArrayList<String> inputArgument = new ArrayList<String>();
    	// get input parameter from BAM/SAM selector
    		for(Iterator<DataRow> it = inData[0].iterator(); it.hasNext(); )
    			inputArgument.add(it.next().getCell(0).toString());

    	// add the input parameter
    	pars.put(" ", StringUtils.join(inputArgument, " "));
    	// return the GUI parameter
		return pars;
	}
	
	
	@Override
	protected BufferedDataTable[] getOutputData(final ExecutionContext exec, String command, final BufferedDataTable[] inData) {
		BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, StringCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, StringCell.TYPE).createSpec()}));
    	
    	DataCell[] c = new DataCell[]{
    			new StringCell(getAbsoluteFilename(SET_OUTPUT_FILE.getStringValue(), true)),
    			new StringCell(command)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
		
        return new BufferedDataTable[]{cont.getTable()};
	}
        
    
    /**
     * Checks, if the file is a annotation file in the GTF / SAF file
     * @param path2AnnotationFile
     * @return
     * @throws InvalidSettingsException
     */
    protected boolean validateAnnotationFile(String path2AnnotationFile) throws InvalidSettingsException {
    	// check for relative path if binary file is valid
    	if(isBinaryValid(getBinaryPath()))
    		path2AnnotationFile = getAbsoluteFilename(path2AnnotationFile, false);
    	
    	// check if path2AnnotationFile exists
    	File f = new File(path2AnnotationFile);
    	if(!(f.isFile() && f.exists()))
    		throw new InvalidSettingsException("Annotation files '" + path2AnnotationFile + "' does not exist.");
    	
    	// TODO: check if GTF or SAF file is valid
    	
    	//	throw new InvalidSettingsException("Folder was found but it seems that it does not contain a valid genome index.");    	
    	// all checks where ok */
    	return true;
    }

	@Override
	protected boolean isParameterEscapingEnabled() {
		return false;
	}
	
	@Override
	protected File getPathToLogOutputFolder() {
		return new File(getAbsoluteFilename(SET_OUTPUT_FILE.getStringValue(), false)).getParentFile();
	}
	
	@Override
	protected File getPathToLockFile() {
		return new File(getAbsoluteFilename(SET_OUTPUT_FILE.getStringValue(), false) + SuccessfulRunChecker.LOCK_ENDING);
	}
}


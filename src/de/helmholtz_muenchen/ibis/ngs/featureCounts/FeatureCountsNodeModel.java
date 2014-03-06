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
    protected static final String CFGKEY_OUTPUT_FILE 		= "OutputFolder";
    protected static final String CFGKEY_ANNOTATION_FILE	= "AnnotationFile";
    protected static final String CFGKEY_ANNOTATION_TYPE	= "AnnotationType";
    protected static final String CFGKEY_ANNOTATION_FEATURE	= "AnnotationFeature";
    protected static final String CFGKEY_THREAD_NUMBER		= "ThreadNumber";
    protected static final String CFGKEY_COUNT_MULTIMAPPED	= "MultimappedFlag";
    
    // initial default values for SettingsModels
    protected static final String DEFAULT_ANNOTATION_TYPE		= "GTF";			// input of GTF file as annotation is default value
    protected static final String ALTERNATIVE_ANNOTATION_TYPE 	= "SAF";			// alternative annotation type
    private static final String DEFAULT_OUTPUT_FOLDER 			= "./output/";		// creates a folder "output" relative to the STAR binary
    private static final String DEFAULT_ANNOTATION_FILE			= "";
    private static final String DEFAULT_ANNOTATION_FEATURE		= "exon";			// default feature which is used for counting
    private static final int DEFAULT_THREAD_NUMBER				= 1;				// default threads to use
    private static final boolean DEAFULT_COUNT_MULTIMAPPED		= false;			// do not count multimapped reads by default

    // name of parameters which are defined in the featureCounts binary
    private final static String NAME_OF_OUTPUT_FILE 		= "-o";		// output file
    private final static String NAME_OF_ANNOTATION_FILE 	= "-a";		// input annotation file
    private final static String NAME_OF_THREAD_NUMBER		= "-T";		// sets the number of threads
    private final static String NAME_OF_FEATURE_TYPE		= "-t";		// sets the name of the feature which is used for counting
    private final static String NAME_OF_COUNT_MULTIMAPPED	= "-M";		// sets if multimapped reads should be counted or not
    private final static String NAME_OF_ANNOTATION_TYPE		= "-F";		// annotation type parameter
    
    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_FEATURE_TYPE			= getSettingsModelString(CFGKEY_ANNOTATION_TYPE);
    private final SettingsModelString SET_OUTPUT_FILE 			= getSettingsModelString(CFGKEY_OUTPUT_FILE);
    private final SettingsModelString SET_ANNOTATION_FILE		= getSettingsModelString(CFGKEY_ANNOTATION_FILE);
    private final SettingsModelString SET_ANNOTATION_TYPE		= getSettingsModelString(CFGKEY_ANNOTATION_TYPE);
    private final SettingsModelInteger SET_THREAD_NUMBER		= getSettingsModelInteger(CFGKEY_THREAD_NUMBER);
    private final SettingsModelBoolean SET_COUNT_MULTIMAPPED	= getSettingsModelBoolean(CFGKEY_COUNT_MULTIMAPPED);
    
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
		if(SET_COUNT_MULTIMAPPED.isEnabled())
			pars.put(NAME_OF_COUNT_MULTIMAPPED, "");
		if(SET_THREAD_NUMBER.getIntValue() > 1)
			pars.put(NAME_OF_THREAD_NUMBER, Integer.toString(SET_THREAD_NUMBER.getIntValue()));
		
		pars.put(NAME_OF_ANNOTATION_TYPE, SET_ANNOTATION_TYPE.getStringValue());
		pars.put(NAME_OF_FEATURE_TYPE, SET_FEATURE_TYPE.getStringValue());
		pars.put(NAME_OF_ANNOTATION_FILE, SET_ANNOTATION_FILE.getStringValue());
		
		/********************* OUTPUT ****************************/
		String outputFolderArgument = getAbsoluteFilename(SET_OUTPUT_FILE.getStringValue(), true);
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
    		path2AnnotationFile = getAbsoluteFilename(path2AnnotationFile, true);
    	
    	// check if genomePath exists
    	File f = new File(path2AnnotationFile);
    	if(!(f.isDirectory() && f.exists()))
    		throw new InvalidSettingsException("Annotation files '" + path2AnnotationFile + "' does not exist.");
    	
    	// TODO: check if GTF or SAF file is valid
    	
    	//	throw new InvalidSettingsException("Folder was found but it seems that it does not contain a valid genome index.");    	
    	// all checks where ok */
    	return true;
    }
}


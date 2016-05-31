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
import org.knime.core.node.util.CheckUtils;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
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
    protected static final String CFGKEY_OUTPUT_FOLDER 				= "OutputFolder";
    protected static final String CFGKEY_ANNOTATION_FILE			= "AnnotationFile";
    protected static final String CFGKEY_ANNOTATION_TYPE			= "AnnotationType";
    protected static final String CFGKEY_ANNOTATION_FEATURE			= "AnnotationFeature";
    protected static final String CFGKEY_THREAD_NUMBER				= "ThreadNumber";
    protected static final String CFGKEY_COUNT_MULTIMAPPED			= "MultimappedFlag";
    protected static final String CFGKEY_COUNT_OVERLAPPING_MULTI	= "MultiOverlapping";
    protected static final String CFGKEY_COUNT_FRAGMENTS			= "Fragments";
    protected static final String CFGKEY_COUNT_CHIMERIC_FRAGMENTS	= "ChimericFragments";
    protected static final String CFGKEY_COUNT_ON_FEATURE_LVL		= "CountOnFeatureLvl";
    protected static final String CFGKEY_GROUP_FEATURE				= "GroupFeature";
    
    // initial default values for SettingsModels
    protected static final String DEFAULT_ANNOTATION_TYPE			= "GTF";			// input of GTF file as annotation is default value
    protected static final String ALTERNATIVE_ANNOTATION_TYPE 		= "SAF";			// alternative annotation type
    protected static final String DEFAULT_ANNOTATION_FILE			= "";
    protected static final String DEFAULT_ANNOTATION_FEATURE		= "exon";			// default feature which is used for counting
    protected static final int DEFAULT_THREAD_NUMBER				= 1;				// default threads to use
    protected static final boolean DEAFULT_COUNT_MULTIMAPPED		= false;			// do not count multimapped reads by default
    protected static final boolean DEAFULT_COUNT_MULTI_OVERLAPING	= false;			// do not count multi overlapping reads by default
    protected static final boolean DEAFULT_COUNT_FRAGMENTS			= false;			// only for paired reads
    protected static final boolean DEAFULT_COUNT_CHIMERIC_FRAGMENTS = false;			// do not count chimeric fragments
    protected static final boolean DEFAULT_COUNT_ON_FEATURE_LVL		= false;			// read summarization will be performed at the feature level (eg. exon level)
    protected static final String DEFAULT_GROUP_FEATURE				= "gene_id";		// group results using id

    
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
    private final static String NAME_OF_COUNT_ON_FEATURE_LVL= "-f";		// enables count on feature level
    private final static String NAME_OF_GROUP_FEATURE		= "-g";		// sets the group id for output file
    				
    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_FEATURE_TYPE				= new SettingsModelString(CFGKEY_ANNOTATION_FEATURE, DEFAULT_ANNOTATION_FEATURE);
    private final SettingsModelString SET_OUTPUT_FOLDER 			= new SettingsModelString(CFGKEY_OUTPUT_FOLDER, "");
    private final SettingsModelString SET_ANNOTATION_FILE			= new SettingsModelString(CFGKEY_ANNOTATION_FILE, DEFAULT_ANNOTATION_FILE);
    private final SettingsModelString SET_ANNOTATION_TYPE			= new SettingsModelString(CFGKEY_ANNOTATION_TYPE, DEFAULT_ANNOTATION_TYPE);
    private final SettingsModelInteger SET_THREAD_NUMBER			= new SettingsModelInteger(CFGKEY_THREAD_NUMBER, DEFAULT_THREAD_NUMBER);
    private final SettingsModelBoolean SET_COUNT_MULTIMAPPED		= new SettingsModelBoolean(CFGKEY_COUNT_MULTIMAPPED, DEAFULT_COUNT_MULTIMAPPED);
    private final SettingsModelBoolean SET_COUNT_OVERLAPPING_MULTI	= new SettingsModelBoolean(CFGKEY_COUNT_OVERLAPPING_MULTI, DEAFULT_COUNT_MULTI_OVERLAPING);
    private final SettingsModelBoolean SET_COUNT_FRAGMENTS			= new SettingsModelBoolean(CFGKEY_COUNT_FRAGMENTS, DEAFULT_COUNT_FRAGMENTS);
    private final SettingsModelBoolean SET_CHIMERIC_FRAGMENTS		= new SettingsModelBoolean(CFGKEY_COUNT_CHIMERIC_FRAGMENTS, DEAFULT_COUNT_CHIMERIC_FRAGMENTS);
    private final SettingsModelBoolean SET_FEATURE_LEVEL			= new SettingsModelBoolean(CFGKEY_COUNT_ON_FEATURE_LVL, DEFAULT_COUNT_ON_FEATURE_LVL);
    private final SettingsModelString SET_GROUP_FEATURE				= new SettingsModelString(CFGKEY_GROUP_FEATURE, DEFAULT_GROUP_FEATURE);
    
    protected final static int MIN_THREADS = 1;
    protected final static int MAX_THREADS = 16;
    
    private String outfile;
    
    // the logger instance
    @SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FeatureCountsNodeModel.class);
       
    /**
     * Constructor for the node model.
     */
    protected FeatureCountsNodeModel() {
        super(1, 1, true, true);
        addSetting(SET_FEATURE_TYPE);
    	addSetting(SET_OUTPUT_FOLDER);
    	addSetting(SET_ANNOTATION_FILE);
    	addSetting(SET_ANNOTATION_TYPE);
    	addSetting(SET_THREAD_NUMBER);
    	addSetting(SET_COUNT_MULTIMAPPED);
    	addSetting(SET_COUNT_OVERLAPPING_MULTI);
    	addSetting(SET_COUNT_FRAGMENTS);
    	addSetting(SET_CHIMERIC_FRAGMENTS);
    	addSetting(SET_FEATURE_LEVEL);
    	addSetting(SET_GROUP_FEATURE);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	super.configure(inSpecs);
    	
    	
    	if(CompatibilityChecker.getIndexCellType(inSpecs[0], "BAMCell")!=0){
    		throw new InvalidSettingsException("Invalid input. No BAMCell in first column of input table. Node requires the FileLoader as predecessor.");
    	}
    	
    	String outfolder_warning = CheckUtils.checkDestinationDirectory(SET_OUTPUT_FOLDER.getStringValue());
		if(outfolder_warning!=null) {
			setWarningMessage(outfolder_warning);
		}
    	
        validateAnnotationFile(SET_ANNOTATION_FILE.getStringValue());

		return new DataTableSpec[]{getDataOutSpec1()};
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
		if(SET_FEATURE_LEVEL.getBooleanValue())
			pars.put(NAME_OF_COUNT_ON_FEATURE_LVL, "");
		if(SET_THREAD_NUMBER.getIntValue() > 1)
			pars.put(NAME_OF_THREAD_NUMBER, Integer.toString(SET_THREAD_NUMBER.getIntValue()));
		
		pars.put(NAME_OF_ANNOTATION_TYPE, SET_ANNOTATION_TYPE.getStringValue());
		pars.put(NAME_OF_FEATURE_TYPE, SET_FEATURE_TYPE.getStringValue());
		pars.put(NAME_OF_ANNOTATION_FILE, SET_ANNOTATION_FILE.getStringValue());
		pars.put(NAME_OF_GROUP_FEATURE, SET_GROUP_FEATURE.getStringValue());
		
		/********************* OUTPUT ****************************/
//		String outputFolderArgument = getAbsoluteFilename(SET_OUTPUT_FOLDER.getStringValue(), false);
//    	File outDir = new File(outputFolderArgument).getParentFile();
//    	// create folder, if not already there
//    	if(!outDir.isDirectory())
//    		outDir.mkdir();
    	
//    	pars.put(NAME_OF_OUTPUT_FILE, outputFolderArgument);

    	
    	/********************** INPUT BAM/SAM ****************************/
    	ArrayList<String> inputArgument = new ArrayList<String>();
    	boolean first = true;
    	String infile;
    	// get input parameter from BAM/SAM selector
    	for(Iterator<DataRow> it = inData[0].iterator(); it.hasNext(); ) {
    		infile = it.next().getCell(0).toString();
    		inputArgument.add(infile);
    		if(first){
				outfile = SET_OUTPUT_FOLDER.getStringValue()+ System.getProperty("file.separator")+ new File(infile).getName();
				outfile = IO.replaceFileExtension(outfile, ".featureCounts");
				first=false;
			}
    	}

    	// add the input parameter
    	pars.put(" ", StringUtils.join(inputArgument, " "));
    	
    	// add the outfile
    	pars.put(NAME_OF_OUTPUT_FILE, outfile);
    	// return the GUI parameter
		return pars;
	}
	
    /**
     * returns the first output specifications of this node
     * @return
     */
    private DataTableSpec getDataOutSpec1() {
    	return new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, StringCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, StringCell.TYPE).createSpec()});
    }
	
	
	@Override
	protected BufferedDataTable[] getOutputData(final ExecutionContext exec, String command, final BufferedDataTable[] inData) {
		BufferedDataContainer cont = exec.createDataContainer(getDataOutSpec1());
		
    	DataCell[] c = new DataCell[]{
    			new StringCell(getAbsoluteFilename(SET_OUTPUT_FOLDER.getStringValue(), true)),
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
	protected File getPathToStderrFile() {
		return new File(outfile + ".out");
	}

	@Override
	protected File getPathToStdoutFile() {
		return new File(outfile + ".err");
	}
		
	@Override
	protected File getPathToLockFile() {
		return new File(outfile + SuccessfulRunChecker.LOCK_ENDING);
	}

}


/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


package de.helmholtz_muenchen.ibis.ngs.mats;

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
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.ngs.featureCounts.FeatureCountsNodeModel;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeModel;

/**
 * This is the model implementation of MATS.
 * 
 *
 * @author Michael Kluge
 */
public class MatsNodeModel extends BinaryWrapperNodeModel {
	
	// name of the output variables
	public static final String OUT_COL1 = "OutputFile";
	public static final String OUT_COL2 = "CallCommand";
    
	// keys for SettingsModels
	protected static final String CFGKEY_OUTPUT_FOLDER 				= "MATS_OutputFolder";
	protected static final String CFGKEY_ANNOTATION_FILE			= "MATS_AnnotationFile";
	protected static final String CFGKEY_READ_LENGTH				= "MATS_ReadLength";
	protected static final String CFGKEY_CUTOFF_DIFFERENCE			= "MATS_CutoffDifference";
	protected static final String CFGKEY_ANALYSIS_TYPE				= "MATS_Analysis_type";
	protected static final String CFGKEY_EXPRESSION_CHANGE			= "MATS_ExpressionChange";
	
	// initial default values for SettingsModels
	protected static final String DEFAULT_OUTPUT_FOLDER 			= "./output/";		// creates a folder "output" relative to the STAR binary
	protected static final String DEFAULT_ANNOTATION_FILE			= "";
	protected static final int DEFAULT_READ_LENGTH					= 100;				// default length of the reads
	protected static final double DEFAULT_CUTOFF_DIFFERENCE			= 0.0001;			// detect also small changes starting at 0.01%
	protected static final boolean DEFAULT_PAIRED_ANALYSIS			= false;			// default is unpaired analysis
	protected static final double DEFAULT_EXPRESSION_CHANGE			= 10000.0;			// exclude AS events when gene expression levels differ more than this foldchange
	
	// name of parameters which are defined in the featureCounts binary
	private final static String NAME_OF_OUTPUT_FOLDER 		= "-o";					// output folder
	private final static String NAME_OF_ANNOTATION_FILE 	= "-gtf";				// input annotation file
	private final static String NAME_OF_READ_LENGTH			= "-len";				// read length
	private final static String NAME_OF_CUTOFF_DIFFERENCE	= "-c";					// cutoff splice difference
	private final static String NAME_OF_PAIRED_ANALYSIS		= "-analysis";			// analysis type
	private final static String NAME_OF_EXPRESSION_CHANGE	= "-expressionChange";	// expression exclusion foldchange 
	
	private final static String NAME_OF_BAM1				= "-b1";	// first condition set 
	private final static String NAME_OF_BAM2				= "-b2";	// second condition set
	private final static String NAME_OF_INSERT_SIZE1		= "-r1";	// insert size of condition set 1
	private final static String NAME_OF_INSERT_SIZE2		= "-r2";	// insert size of condition set 2
	private final static String NAME_OF_STANDARD_DEV1		= "-sd1";	// standard deviation of condition set 1 
	private final static String NAME_OF_STANDARD_DEV2		= "-sd2";	// standard deviation of condition set 2  
	
	
	// definition of SettingsModel (all prefixed with SET)
	private final SettingsModelString SET_OUTPUT_FOLDER			= new SettingsModelString(CFGKEY_OUTPUT_FOLDER, DEFAULT_OUTPUT_FOLDER);
	private final SettingsModelString SET_ANNOTATION_FILE		= new SettingsModelString(CFGKEY_ANNOTATION_FILE, DEFAULT_ANNOTATION_FILE);
	private final SettingsModelInteger SET_READ_LENGTH			= new SettingsModelInteger(CFGKEY_READ_LENGTH, DEFAULT_READ_LENGTH);	
	private final SettingsModelDouble SET_CUTOFF_DIFFERENCE		= new SettingsModelDouble(CFGKEY_CUTOFF_DIFFERENCE, DEFAULT_CUTOFF_DIFFERENCE);
	private final SettingsModelBoolean SET_ANALYSIS_TYPE		= new SettingsModelBoolean(CFGKEY_ANALYSIS_TYPE, DEFAULT_PAIRED_ANALYSIS);
	private final SettingsModelDouble SET_EXPRESSION_CHANGE		= new SettingsModelDouble(CFGKEY_EXPRESSION_CHANGE, DEFAULT_EXPRESSION_CHANGE);
	
	protected final static double VALID_FOLDCHANGE_CUTOFF = 1.0;
	protected final static double MIN_CHANGE = 0.0;
	protected final static double MAX_CHANGE = 1.0;
	
	// the logger instance
	@SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FeatureCountsNodeModel.class);


    /**
     * Constructor for the node model.
     */
    protected MatsNodeModel() {
    	super(2, 1, true, true);
    	addSetting(SET_OUTPUT_FOLDER);
    	addSetting(SET_ANNOTATION_FILE);
    	addSetting(SET_READ_LENGTH);
    	addSetting(SET_CUTOFF_DIFFERENCE);
    	addSetting(SET_ANALYSIS_TYPE);
    	addSetting(SET_EXPRESSION_CHANGE);
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
		if(SET_ANALYSIS_TYPE.getBooleanValue())
			pars.put(NAME_OF_PAIRED_ANALYSIS, "P");
		else
			pars.put(NAME_OF_PAIRED_ANALYSIS, "U");
		if(SET_READ_LENGTH.getIntValue() > 1)
			pars.put(NAME_OF_READ_LENGTH, Integer.toString(SET_READ_LENGTH.getIntValue()));
		
		pars.put(NAME_OF_ANNOTATION_FILE, SET_ANNOTATION_FILE.getStringValue());
		pars.put(NAME_OF_CUTOFF_DIFFERENCE, Double.toString(Math.min(SET_CUTOFF_DIFFERENCE.getDoubleValue(), VALID_FOLDCHANGE_CUTOFF)));
		pars.put(NAME_OF_EXPRESSION_CHANGE, Double.toString(SET_EXPRESSION_CHANGE.getDoubleValue()));
		
		
		/********************* OUTPUT ****************************/
		String outputFolderArgument = getAbsoluteFilename(SET_OUTPUT_FOLDER.getStringValue(), true);
    	File outDir = new File(outputFolderArgument).getParentFile();
    	// create folder, if not already there
    	if(!outDir.isDirectory())
    		outDir.mkdir();
    	
    	pars.put(NAME_OF_OUTPUT_FOLDER, outputFolderArgument);

    	
    	/********************** INPUT BAM SETTINGS OF CONDITION 1 ****************************/
    	ArrayList<String> bamArgument = new ArrayList<String>();
    	ArrayList<String> insertArgument = new ArrayList<String>();
    	ArrayList<String> standardDevArgument = new ArrayList<String>();
    	
    	// get input parameter from BAM/SAM selector
		for(Iterator<DataRow> it = inData[0].iterator(); it.hasNext(); ) {
			DataRow row = it.next();
			bamArgument.add(row.getCell(0).toString());
			insertArgument.add(row.getCell(1).toString());
			standardDevArgument.add(row.getCell(2).toString());
		}
    	
    	// add the input parameter
		pars.put(NAME_OF_BAM1, StringUtils.join(bamArgument, ","));
		pars.put(NAME_OF_INSERT_SIZE1, StringUtils.join(insertArgument, ","));
		pars.put(NAME_OF_STANDARD_DEV1, StringUtils.join(standardDevArgument, ","));
		
		/********************** INPUT BAM SETTINGS OF CONDITION 1 ****************************/
		bamArgument.clear();
		insertArgument.clear();
		standardDevArgument.clear();
		
		// get input parameter from BAM/SAM selector
		for(Iterator<DataRow> it = inData[1].iterator(); it.hasNext(); ) {
			DataRow row = it.next();
			bamArgument.add(row.getCell(0).toString());
			insertArgument.add(row.getCell(1).toString());
			standardDevArgument.add(row.getCell(2).toString());
		}
		
    	// add the input parameter
		pars.put(NAME_OF_BAM2, StringUtils.join(bamArgument, ","));
		pars.put(NAME_OF_INSERT_SIZE2, StringUtils.join(insertArgument, ","));
		pars.put(NAME_OF_STANDARD_DEV2, StringUtils.join(standardDevArgument, ","));
		
    	// return the GUI parameter
		return pars;
	}
	
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
		// validate input data
		validateInputTables(inData);
		return super.execute(inData, exec);
	}
	
	/**
	 * Validates the input tables
	 * @param inData
	 * @return
	 * @throws InvalidSettingsException 
	 */
    protected boolean validateInputTables(final BufferedDataTable[] inData) throws InvalidSettingsException {
    	if(inData.length == 2) {
    		for(int i = 0; i < 2; i++) {
    			for(Iterator<DataRow> it = inData[i].iterator(); it.hasNext(); ) {
    				// check, if first one is file
    				String fname = it.next().getCell(0).toString();
    				if(!new File(fname).exists()) {
    					throw new InvalidSettingsException("BAM file '" + fname + "' does not exist."); 
    				}
    			}
    		}
    	}
    	
    	return true;
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
    	
    	// TODO: check if GTF file is valid	
    	// all checks where ok */
    	return true;
    }
    
	@Override
	protected File getPathToStderrFile() {
		return new File(getAbsoluteFilename(SET_OUTPUT_FOLDER.getStringValue(), true) + File.separator + "log.out");
	}

	@Override
	protected File getPathToStdoutFile() {
		return new File(getAbsoluteFilename(SET_OUTPUT_FOLDER.getStringValue(), true) + File.separator + "log.err");
	}
		
	@Override
	protected File getPathToLockFile() {
		return new File(getAbsoluteFilename(SET_OUTPUT_FOLDER.getStringValue(), true) + File.separator + "LOCK" + SuccessfulRunChecker.LOCK_ENDING);
	}
}

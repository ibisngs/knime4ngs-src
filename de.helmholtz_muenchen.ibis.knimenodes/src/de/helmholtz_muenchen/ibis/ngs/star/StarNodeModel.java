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
package de.helmholtz_muenchen.ibis.ngs.star;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.SAMCell;

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
	public static final String OUT_COL = "Output";
	
    // keys for SettingsModels
    protected static final String CFGKEY_RUN_MODE 		= "RunMode";
    protected static final String CFGKEY_OUTPUT_FOLDER 	= "OutputFolder";
    protected static final String CFGKEY_GENOME_FOLDER 	= "GenomeFolder";
    protected static final String CFGKEY_OPTIONAL_PARA 	= "OPTIONAL_PARA";

    // initial default values for SettingsModels    
    protected static final String DEFAULT_RUN_MODE 		= "alignReads";			// align read mode is default
    protected static final String ALTERNATIVE_RUN_MODE 	= "genomeGenerate";		// alternative run mode
    protected static final String DEFAULT_OUTPUT_FOLDER 	= "";			// creates a folder "output" relative to the STAR binary
    protected static final String DEFAULT_GENOME_FOLDER 	= "";		// searches for the genome index in "GenomeDir" relative to the STAR binary

    // name of parameters which are defined in the STAR binary
    private final static String NAME_OF_OUTPUT_PREFIX_PARAM 	= "--outFileNamePrefix";	// parameter in STAR which allows the user to set the output folder (relative or absolute)
    private final static String NAME_OF_OUTPUT_GENOMEDIR_PARAM 	= "--genomeDir";			// output parameter in case of generateGenome run or input in other case
    private final static String NAME_OF_RUN_MODE_PARAM 			= "--runMode";				// parameter in STAR which sets the runMode
    private final static String NAME_OF_FASTAQ_PARAM 			= "--readFilesIn";			// string(s): paths to files that contain input read1 (and, if needed,  read2)
    private final static String NAME_OF_FASTA_FILES_PARAM		= "--genomeFastaFiles";		// fasta files with genomic sequence for genome files generation
    
    private final static String NAME_OF_GENOME_PARAMETER_FILE	= "genomeParameters.txt"; 	// name of settings file from a indexed genome 
	
    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_RUN_MODE					= new SettingsModelString(CFGKEY_RUN_MODE, DEFAULT_RUN_MODE);
    private final SettingsModelString SET_OUTPUT_FOLDER				= new SettingsModelString(CFGKEY_OUTPUT_FOLDER, DEFAULT_OUTPUT_FOLDER);
    private final SettingsModelString SET_GENOME_FOLDER				= new SettingsModelString(CFGKEY_GENOME_FOLDER, DEFAULT_GENOME_FOLDER);
    private final SettingsModelOptionalString SET_OPTIONAL_PARA		= new SettingsModelOptionalString(CFGKEY_OPTIONAL_PARA, "",false);
    private String OUTFILE 											= "";
//    private String OUTFILE_TABLE									= "";
    
    // the logger instance
    @SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(StarNodeModel.class);
	private static String readType = "";
       
    /**
     * Constructor for the node model.
     */
    protected StarNodeModel() {
        super(1, 1, true, true);
        addSetting(SET_RUN_MODE);
    	addSetting(SET_OUTPUT_FOLDER);
    	addSetting(SET_GENOME_FOLDER);
    	addSetting(SET_OPTIONAL_PARA);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	super.configure(inSpecs);

    	CompatibilityChecker CC = new CompatibilityChecker();
    	
    	if(isAlignRunMode()){
    		readType = CC.getReadType(inSpecs, 0);
        	if(CC.getWarningStatus()){
        		setWarningMessage(CC.getWarningMessages());
        	}
        	validateGenomeIndex(SET_GENOME_FOLDER.getStringValue());
    	} else {
    		if(!CompatibilityChecker.checkInputCellType(inSpecs[0], "FastACell")) {
    			throw new InvalidSettingsException("Incompatible input: In 'genomeGenerate' mode the node expects a FastA file as input.");

    		}
    	}
    	
		return new DataTableSpec[]{getDataOutSpec1()};
    }
    
    

	@Override
	protected LinkedHashMap<String, String> getGUIParameters(final BufferedDataTable[] inData) {
		LinkedHashMap<String, String> pars = new LinkedHashMap<String, String>();

		/********************* RUN MODE ***************************/
		pars.put(NAME_OF_RUN_MODE_PARAM, SET_RUN_MODE.getStringValue());
		
		   	
    	/********************** INPUT ****************************/
    	String inputParameter = (isAlignRunMode() ? NAME_OF_FASTAQ_PARAM : NAME_OF_FASTA_FILES_PARAM);
    	ArrayList<String> inputArgument = new ArrayList<String>();
    	
    	// get input parameter from run aligner
    	if(isAlignRunMode()) {	
	    	String path2readFile1 = inData[0].iterator().next().getCell(0).toString();
	    	String path2readFile2 = "";
	    	
	    	// add first input file
	    	inputArgument.add(getAbsoluteFilename(path2readFile1, false));

	    	if(readType.equals("paired-end")) {
	    		path2readFile2 = inData[0].iterator().next().getCell(1).toString();
		    	// add second input file, if paired mode was selected and both files are different
		    	if(!path2readFile1.equals(path2readFile2) && path2readFile2.length() > 0)
		    		inputArgument.add(getAbsoluteFilename(path2readFile2, false));
	    	}
	    	
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
    	
    	//add optional parameters
    	if(SET_OPTIONAL_PARA.isActive()){
    		pars.put(SET_OPTIONAL_PARA.getStringValue(), "");
    	}
    	
    	
    	
		/********************* OUTPUT ****************************/
    	// check, which kind of output parameter must be set.
    	String infile = inData[0].iterator().next().getCell(0).toString();
    	String outputFolderParameter = (isAlignRunMode() ? NAME_OF_OUTPUT_PREFIX_PARAM : NAME_OF_OUTPUT_GENOMEDIR_PARAM);
		String outputFolderArgument = SET_OUTPUT_FOLDER.getStringValue();
		if(CompatibilityChecker.inputFileNotOk(outputFolderArgument, false)) {
			outputFolderArgument = new File(infile).getParent();
			if(!isAlignRunMode()) {
				outputFolderArgument += File.separator + "STAR_genome";
			}
		}
		
		if(!outputFolderArgument.endsWith(File.separator)) {
			outputFolderArgument += File.separator;
		}
		
    	File outDir = new File(outputFolderArgument);
    	// create folder, if not already there
    	if(!outDir.isDirectory())
    		outDir.mkdirs();
    	
    	if(isAlignRunMode()){
    		String outfile = IO.replaceFileExtension(infile,".");
    		outfile = IO.getFileName(outfile);
    		outputFolderArgument+=outfile;	
    		
        	pars.put(outputFolderParameter, outputFolderArgument);
        	OUTFILE = outputFolderArgument+"Aligned.out.sam";
    	}else{
    		pars.put(NAME_OF_OUTPUT_PREFIX_PARAM, outputFolderArgument);
        	pars.put(outputFolderParameter, outputFolderArgument);
        	OUTFILE = outputFolderArgument;
    	}
	
    	// return the GUI parameter
		return pars;
	}
	
    /**
     * returns the first output specifications of this node
     * @return
     */
    private DataTableSpec getDataOutSpec1() {
    	
    	if(isAlignRunMode()) {
        	return new DataTableSpec(
        			new DataColumnSpec[]{
        					new DataColumnSpecCreator(OUT_COL, SAMCell.TYPE).createSpec()});
    	}else{
        	return new DataTableSpec(
        			new DataColumnSpec[]{
        					new DataColumnSpecCreator(OUT_COL, FileCell.TYPE).createSpec()});
    	}
    }
	
	@Override
	protected BufferedDataTable[] getOutputData(final ExecutionContext exec, String command, final BufferedDataTable[] inData) {
		BufferedDataContainer cont = exec.createDataContainer(getDataOutSpec1());

		String tmpOut;
		if(isAlignRunMode()) {
			tmpOut = OUTFILE.replaceFirst("Aligned.out.sam", "_STARtmp");	
		} else {
			tmpOut = OUTFILE + "_STARtmp";
		}
		
		if(tmpOut.contains("_STARtmp")) {
			File f = new File(tmpOut);
			if(f.exists()){
				try{
					FileUtils.deleteDirectory(f);
				}catch(Exception e){
					setWarningMessage("Failed to delete tmp dir.");
				}
			}
		}

    	DataCell[] c = new FileCell[]{ FileCellFactory.create(OUTFILE)};
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
    

	@Override
	protected File getPathToLockFile() {
		if(!isAlignRunMode()) return new File(OUTFILE+"Genome"+ SuccessfulRunChecker.LOCK_ENDING);
		return new File(OUTFILE + SuccessfulRunChecker.LOCK_ENDING);
	}

	@Override
	protected File getPathToStderrFile() {
		if(!isAlignRunMode()) return new File(OUTFILE+"Genome.stdErr");
		return new File(OUTFILE + ".stdErr");
	}

	@Override
	protected File getPathToStdoutFile() {
		if(!isAlignRunMode()) return new File(OUTFILE+"Genome.stdOut");
		return new File(OUTFILE + ".stdOut");
 	}
}


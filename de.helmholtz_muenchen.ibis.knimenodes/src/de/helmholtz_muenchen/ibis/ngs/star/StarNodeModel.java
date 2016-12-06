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
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
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
     
    public static final String CFGKEY_RUN_MODE 			= "RunMode";
    public static final String CFGKEY_TWOPASS			= "two_pass_mode";
    public static final String CFGKEY_OUTPUT_FOLDER 	= "OutputFolder";
    public static final String CFGKEY_GENOME_FOLDER 	= "GenomeFolder";
    public static final String CFGKEY_GTF_FILE			= "gtf_file";
    public static final String CFGKEY_OVERHANG			= "overhang";
    public static final String CFGKEY_THREADS			= "threads";
    public static final String CFGKEY_OPTIONAL_PARA 	= "opt_parameter";
   
    public static final String DEFAULT_RUN_MODE 		= "alignReads";
    public static final String ALTERNATIVE_RUN_MODE 	= "genomeGenerate";

    private final SettingsModelString m_runMode					= new SettingsModelString(CFGKEY_RUN_MODE, DEFAULT_RUN_MODE);
    private final SettingsModelBoolean m_two_pass				= new SettingsModelBoolean(CFGKEY_TWOPASS, true);
    private final SettingsModelString m_outfolder				= new SettingsModelString(CFGKEY_OUTPUT_FOLDER, "");
    private final SettingsModelString m_genome_folder			= new SettingsModelString(CFGKEY_GENOME_FOLDER, "");
    private final SettingsModelString m_gtf_file				= new SettingsModelString(CFGKEY_GTF_FILE, "");
    private final SettingsModelIntegerBounded m_overhang		= new SettingsModelIntegerBounded(CFGKEY_OVERHANG, 100, 1, Integer.MAX_VALUE);
    private final SettingsModelIntegerBounded m_threads			= new SettingsModelIntegerBounded(CFGKEY_THREADS, 4, 1, Integer.MAX_VALUE);
    private final SettingsModelOptionalString m_opt_parameter	= new SettingsModelOptionalString(CFGKEY_OPTIONAL_PARA, "",false);
    
	private String readType, genome_folder, out_folder, gtf_file, OUTFILE;
	
	public static final String OUT_COL = "Output";
       
    /**
     * Constructor for the node model.
     */
    protected StarNodeModel() {
        super(1, 1, true, true);
        addSetting(m_runMode);
        addSetting(m_two_pass);
    	addSetting(m_outfolder);
    	addSetting(m_genome_folder);
    	addSetting(m_threads);
    	addSetting(m_gtf_file);
    	addSetting(m_opt_parameter);
    	
    	addPrefPageSetting(SET_BINARY_PATH,IBISKNIMENodesPlugin.STAR);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	super.updatePrefs();
    	
    	if(CompatibilityChecker.inputFileNotOk(getBinaryPath(), false)) {
			throw new InvalidSettingsException("Set path to featureCounts binary!");
		}
    	
    	genome_folder = IO.processFilePath(m_genome_folder.getStringValue());
    	out_folder = IO.processFilePath(m_outfolder.getStringValue());
    	gtf_file = IO.processFilePath(m_gtf_file.getStringValue());
    	
    	CompatibilityChecker CC = new CompatibilityChecker();
    	if(isAlignRunMode()){
    		readType = CC.getReadType(inSpecs, 0);
        	if(CC.getWarningStatus()){
        		setWarningMessage(CC.getWarningMessages());
        	}
        	validateGenomeIndex(genome_folder);
    	} else {
    		if(!CompatibilityChecker.checkInputCellType(inSpecs[0], "FastACell")) {
    			throw new InvalidSettingsException("Incompatible input: In 'genomeGenerate' mode the node expects a FastA file as input.");
    		}
    		if(CompatibilityChecker.inputFileNotOk(gtf_file,true)) {
    			setWarningMessage("Using a GTF file is strongly recommended for genomeGenerate!");
    		}
    	}
		return new DataTableSpec[]{getDataOutSpec1()};
    }

	@Override
	protected LinkedHashMap<String, String> getGUIParameters(final BufferedDataTable[] inData) {
		
		LinkedHashMap<String, String> pars = new LinkedHashMap<String, String>();
		pars.put("--runMode", m_runMode.getStringValue());
		  	
		//input files
    	String inputParameter = (isAlignRunMode() ? "--readFilesIn" : "--genomeFastaFiles");
    	ArrayList<String> inputArgument = new ArrayList<String>();
    	
    	boolean is_zipped = false;
    	
    	if(isAlignRunMode()) {	
	    	String path2readFile1 = inData[0].iterator().next().getCell(0).toString();
	    	is_zipped = path2readFile1.endsWith(".gz");
	    	inputArgument.add(path2readFile1);

	    	if(readType.equals("paired-end")) {
	    		String path2readFile2 = inData[0].iterator().next().getCell(1).toString();
		    	if(!path2readFile1.equals(path2readFile2) && path2readFile2.length() > 0 && path2readFile2.endsWith(".gz") == is_zipped) {
		    		inputArgument.add(path2readFile2);
		    	}
	    	}
	    	
	    	if(is_zipped) {
	    		pars.put("--readFilesCommand", "zcat");
	    	}
	    	pars.put("--genomeDir", genome_folder);
	    	
	    	if(m_two_pass.getBooleanValue()) {
	    		pars.put("--twopassMode", "Basic");
	    	}
    	} else {
    		for(DataRow dr : inData[0]) {
    			inputArgument.add(dr.getCell(0).toString());
    		}
    	}
    	pars.put(inputParameter, StringUtils.join(inputArgument, " "));
    	
    	//output parameters
    	String infile = inData[0].iterator().next().getCell(0).toString();
		if(CompatibilityChecker.inputFileNotOk(out_folder, false)) {
			out_folder = new File(infile).getParent();
			if(!isAlignRunMode()) {
				out_folder += File.separator + "STAR_genome";
			}
		}
		
		if(!out_folder.endsWith(File.separator)) {
			out_folder += File.separator;
		}
		
    	File outDir = new File(out_folder);
    	if(!outDir.isDirectory()) {
    		outDir.mkdirs();
    	}
    	
    	OUTFILE = out_folder;
    	if(isAlignRunMode()){
    		String outfile = IO.replaceFileExtension(infile,".");
    		outfile = IO.getFileName(outfile);
    		out_folder+=outfile;	
        	OUTFILE = out_folder+"Aligned.out.sam";
    	}
    	
    	if(!isAlignRunMode()) {
    		pars.put("--genomeDir", out_folder);
    	}
    	
    	//required for log files
    	pars.put("--outFileNamePrefix", out_folder);
    	
    	//add further parameters
    	
    	pars.put("--runThreadN", m_threads.getIntValue()+"");
    	
    	if(!CompatibilityChecker.inputFileNotOk(gtf_file,true)) {
    		pars.put("--sjdbGTFfile", gtf_file);
    		pars.put("--sjdbOverhang", m_overhang.getIntValue()+"");
    	}
    	
    	if(m_opt_parameter.isActive()){
    		pars.put(m_opt_parameter.getStringValue(), "");
    	}
    	
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
    	return DEFAULT_RUN_MODE.equals(m_runMode.getStringValue());
    }
    
    
    /**
     * Checks, if the path contains a genome index for STAR
     * @param genomePath
     * @return
     * @throws InvalidSettingsException
     */
    protected boolean validateGenomeIndex(String genomePath) throws InvalidSettingsException {
    	// check if genomePath exists
    	File f = new File(genomePath);
    	if(!(f.isDirectory() && f.exists()))
    		throw new InvalidSettingsException("Genome path '" + genomePath + "' does not exist.");
    	
    	// check if parameter file is there
    	f = new File(genomePath + File.separator + "genomeParameters.txt");
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

	@Override
	protected String getOutfile() {
		return OUTFILE;
	}
}


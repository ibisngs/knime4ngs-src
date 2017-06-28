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
package de.helmholtz_muenchen.ibis.ngs.trimgalore;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;

import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;


/**
 * This is the model implementation of TrimGalore.
 * Trim Galore! is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing).
 *
 * @author Paul Hager
 */
public class TrimGaloreNodeModel extends HTExecutorNodeModel {
    
    // the logger instance
    private static final NodeLogger LOGGER = NodeLogger
            .getLogger(TrimGaloreNodeModel.class);
        
    // Node specific params
    public static final String CFGKEY_FASTQC = "Path2FastQC";
    public static final String CFGKEY_TRIMG = "Path2TrimGalore";
    public static final String CFGKEY_CUTADAPT = "Path2Cutadapt";
    public static final String CFGKEY_OUTFOLDER_FASTQC = "Path2OutFolderFastQC";
    public static final String CFGKEY_OUTFOLDER_TRIMG = "Path2OutFolderTrimGalore";
    public static final String CFGKEY_FASTQC_ENABLE = "EnableFastQC";
	public static final String CFGKEY_THREADS_FASTQC = "NumberThreads";
	public static final String CFGKEY_QUALITY = "QualityTrimThreshold";
	public static final String CFGKEY_ADAPTER = "AdapterSequence";
	public static final String CFGKEY_ADAPTER2 = "AdapterSequence2";
	public static final String CFGKEY_PRESET_ADAPTER = "PresetAdapter";
	public static final String CFGKEY_MAX_LENGTH = "MaxLength";
	public static final String CFGKEY_STRINGENCY = "Stringency";
	public static final String CFGKEY_ERROR_RATE = "ErrorRate";
	public static final String CFGKEY_GZIP = "GzipCompress";
	public static final String CFGKEY_LENGTH = "Length";
	public static final String CFGKEY_ADDITIONAL_OPTIONS = "AdditionalOptions";
	public static final String CFGKEY_FASTQC_ADDITIONAL_OPTIONS = "FastqcAdditionalOptions";
	
	protected final static int DEFAULT_FASTQC_THREADS = 4;
	protected final static int DEFAULT_QUALITY = 20;
	protected final static int DEFAULT_STRINGENCY = 1;
	protected final static double DEFAULT_ERROR_RATE = 0.1;
	protected final static int DEFAULT_LENGTH = 20;
	
	
	private final SettingsModelString m_fastqc = 
			new SettingsModelString(CFGKEY_FASTQC,"");
	
	private final SettingsModelString m_trimg =
			new SettingsModelString(CFGKEY_TRIMG, "");
	
	private final SettingsModelString m_cutadapt =
			new SettingsModelString(CFGKEY_CUTADAPT, "");
	
	private final SettingsModelString m_outfolder_fastqc = 
			new SettingsModelString(CFGKEY_OUTFOLDER_FASTQC, "");
	
	private final SettingsModelString m_outfolder_trimg = 
			new SettingsModelString(CFGKEY_OUTFOLDER_TRIMG, "");
	
	private final SettingsModelBoolean m_fastqc_enabled =
			new SettingsModelBoolean(CFGKEY_FASTQC_ENABLE, false);
	
	private final SettingsModelIntegerBounded m_threads = 
			new SettingsModelIntegerBounded(CFGKEY_THREADS_FASTQC, DEFAULT_FASTQC_THREADS, 1, Integer.MAX_VALUE);
	
	private final SettingsModelIntegerBounded m_quality = 
			new SettingsModelIntegerBounded(CFGKEY_QUALITY, DEFAULT_QUALITY, 0, Integer.MAX_VALUE);
	
	private final SettingsModelString m_adapter = 
			new SettingsModelString(CFGKEY_ADAPTER, "");
	
	private final SettingsModelString m_adapter2 = 
			new SettingsModelString(CFGKEY_ADAPTER2,"");
	
	private final SettingsModelString m_preset_adapter = 
			new SettingsModelString(CFGKEY_PRESET_ADAPTER,"");
	
	//private final SettingsModelIntegerBounded m_max_length = 
	//		new SettingsModelIntegerBounded(CFGKEY_MAX_LENGTH, 0, 0, Integer.MAX_VALUE);
	
	private final SettingsModelIntegerBounded m_stringency = 
			new SettingsModelIntegerBounded(CFGKEY_STRINGENCY, DEFAULT_STRINGENCY, 1, Integer.MAX_VALUE);
	
	private final SettingsModelDoubleBounded m_error_rate = 
			new SettingsModelDoubleBounded(CFGKEY_ERROR_RATE, DEFAULT_ERROR_RATE, 0, 1);
	
	private final SettingsModelBoolean m_gzip = 
			new SettingsModelBoolean(CFGKEY_GZIP, true);
	
	private final SettingsModelIntegerBounded m_length = 
			new SettingsModelIntegerBounded(CFGKEY_LENGTH, DEFAULT_LENGTH, 0, Integer.MAX_VALUE);
	
	private final SettingsModelString m_additional_options =
			new SettingsModelString(CFGKEY_ADDITIONAL_OPTIONS, "");
	
	private final SettingsModelString m_fastqc_additional_options =
			new SettingsModelString(CFGKEY_FASTQC_ADDITIONAL_OPTIONS, "");
	
	//The Output Col Names
	public static final String OUT_TRIMREAD1 = "Path2TrimmedReadsFile1";
	public static final String OUT_TRIMREAD2 = "Path2TrimmedReadsFile2";
	
	private int count = 2;
	
	//ReadType: paired-end or single-end
	private static String readType = "";
    

    /**
     * Constructor for the node model.
     */
    protected TrimGaloreNodeModel() {
    	
    	
    	super(1, 1);
        
        addSetting(m_fastqc);
        addSetting(m_trimg);
        addSetting(m_cutadapt);
    	addSetting(m_outfolder_fastqc);
    	addSetting(m_outfolder_trimg);
    	addSetting(m_fastqc_enabled);
    	addSetting(m_threads);
    	addSetting(m_quality);
    	addSetting(m_adapter);
    	addSetting(m_adapter2);
    	addSetting(m_preset_adapter);
    	//addSetting(m_max_length);
    	addSetting(m_stringency);
    	addSetting(m_error_rate);
    	addSetting(m_gzip);
    	addSetting(m_length);
    	addSetting(m_additional_options);
    	addSetting(m_fastqc_additional_options);
    	
    	addPrefPageSetting(m_fastqc, IBISKNIMENodesPlugin.FASTQC);
    	addPrefPageSetting(m_trimg, IBISKNIMENodesPlugin.TRIMG);
    	addPrefPageSetting(m_cutadapt, IBISKNIMENodesPlugin.CUTADAPT);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	String cutadaptPath = m_cutadapt.getStringValue();
    	String fastqcPath = m_fastqc.getStringValue();
    	String outFolderTrimG = m_outfolder_trimg.getStringValue();
    	int quality = m_quality.getIntValue();
    	String adapter = m_adapter.getStringValue();
    	String adapter2 = m_adapter2.getStringValue();
    	String presetAdapter = m_preset_adapter.getStringValue();
    	int stringency = m_stringency.getIntValue(); 
    	double errorRate = m_error_rate.getDoubleValue();
    	int length = m_length.getIntValue();
    	String outfolderFastQC = m_outfolder_fastqc.getStringValue();
    	int fastqcThreads = m_threads.getIntValue();
    	String fastqcAddOpts = m_fastqc_additional_options.getStringValue();
    	String trimgAddOpts = m_additional_options.getStringValue();
    	
    	// Beautify outfolder paths
    	outFolderTrimG = outFolderTrimG.trim();
    	if(!outFolderTrimG.equals("") && !outFolderTrimG.endsWith(System.getProperty("file.separator"))){
    		outFolderTrimG += System.getProperty("file.separator");
    	}
    	
    	outfolderFastQC = outfolderFastQC.trim();
    	if(!outfolderFastQC.equals("") && !outfolderFastQC.endsWith(System.getProperty("file.separator"))){
    		outfolderFastQC += System.getProperty("file.separator");
    	}
    	
    	//Check input table integrity
    	CompatibilityChecker.inDataCheck(inData);
    	
    	String inFile1 = inData[0].iterator().next().getCell(0).toString();
    	String inFile2 = "";
    	
    	if(CompatibilityChecker.inputFileNotOk(inFile1)) {
    		throw new InvalidSettingsException("Input FastQ file does not exist or is invalid!");
    	}
    	
    	String path_variable = "PATH=";
    	
    	String outFile = this.createOutputName(inFile1);
    	String outFile2 = "";
    	
    	
    	ArrayList<String> cmd = new ArrayList<String>();
    	cmd.add(IO.processFilePath(m_trimg.getStringValue()));
    	
    	if(cutadaptPath.equals("") || Files.notExists(Paths.get(cutadaptPath))) {
    		setWarningMessage("Cutadapt PATH was not specified!");
    	} else {
    		cmd.add("--path_to_cutadapt="+IO.processFilePath(cutadaptPath));
    	}
    	
    	if(fastqcPath.equals("")|| Files.notExists(Paths.get(fastqcPath))) {
    		setWarningMessage("FastQC PATH was not specified!");
    	} else {
    		//remove 'fastqc' at end of path
    		int pos = fastqcPath.lastIndexOf(System.getProperty("file.separator"));
    		path_variable += fastqcPath.substring(0,pos)+":";
    	}
    	
    	path_variable += System.getenv("PATH");
    	
    	if(!outFolderTrimG.equals("")){
    		cmd.add("-o="+IO.processFilePath(outFolderTrimG));
    	} else {
    		String inFolder = inFile1.substring(0, inFile1.lastIndexOf(System.getProperty("file.separator"))+1);
    		cmd.add("-o="+inFolder);
    	}
    	
    	File lockFile = new File(outFile.substring(0,outFile.lastIndexOf(".")) + ".TrimG" + SuccessfulRunChecker.LOCK_ENDING);
    	
    	if(quality != 20){
    		cmd.add("-q="+quality);
    	}
    	
    	if(!adapter.equals("")){
    		cmd.add("-a="+adapter);
    	}
    	
    	if(!adapter2.equals("")){
    		cmd.add("-a2="+adapter2);
    	}
    	
    	if(!presetAdapter.equals("") && adapter.equals("")){
    		if(presetAdapter.equals("illumina")){
    			cmd.add("--illumina");
    		} else if(presetAdapter.equals("nextera")){
    			cmd.add("--nextera");
    		} else if (presetAdapter.equals("small_rna")){
    			cmd.add("--small_rna");
    		}
    	}
    	
    	/*
    	if(m_max_length.getIntValue() != 0){
    		cmd.add("--max_length="+m_max_length.getIntValue());
    	}
    	*/
    	
    	if(stringency != 1){
    		cmd.add("--stringency="+stringency);
    	}
    	
    	if(errorRate != 0.1){
    		cmd.add("-e="+errorRate);
    	}
    	
    	if(m_gzip.getBooleanValue()){
    		cmd.add("--gzip");
    	} else {
    		cmd.add("--dont_gzip");
    	}
    	
    	if(length != 20){
    		cmd.add("--length="+length);
    	}
    	
    	if(m_fastqc_enabled.getBooleanValue()){
    		String fastqcString = "";
    		if(fastqcThreads != 1 || !outfolderFastQC.equals("")){
    			fastqcString += "--fastqc_args=";
    			if(fastqcThreads != 1){
    				fastqcString += "-t="+fastqcThreads+" ";
    			}
    			if(!outfolderFastQC.equals("")){
    				fastqcString += "-o="+IO.processFilePath(outfolderFastQC)+" ";
    			}
    			if(!fastqcAddOpts.equals("")){
    				setWarningMessage("NOTE! All additional FastQC options MUST use equals operator to bind arguments to options. i.e. --k=5");
    				fastqcString += fastqcAddOpts;	
    			}
    			fastqcString = fastqcString.trim();
    			fastqcString += "";
    		} else {
    			fastqcString += "--fastqc";
    		}
    		cmd.add(fastqcString);
    	}
    	
    	if(!trimgAddOpts.equals("")){
    		setWarningMessage("NOTE! All additional TrimGalore options MUST use equals operator to bind arguments to options. i.e. --length_1=30");
    		String[] addOpts = trimgAddOpts.split(" ");
    		for(String opt : addOpts){
    			cmd.add(opt);
    		}
    	}
    	
    	
    	// Add paired-end options if necessary
    	if(readType.equals("paired-end")){
    		cmd.add("--paired");
    		inFile2 = inData[0].iterator().next().getCell(1).toString();
    		
    		if(CompatibilityChecker.inputFileNotOk(inFile2)) {
        		throw new InvalidSettingsException("Input FastQ file does not exist or is invalid!");
        	}
    		
    		outFile2 = createOutputName(inFile2);
    		
    		cmd.add(inFile2);
    	}
    	
    	cmd.add(inFile1);
    	
    	// Parse options into command string array - ensures spaces in path will be correctly interpreted
    	String[] command = new String[cmd.size()];
    	for(int indx=0; indx < cmd.size(); indx++){
    		command[indx] = cmd.get(indx);
    	}
    	
    	LOGGER.info("CMD: "+command.toString());
    	
    	LOGGER.info("Running Trim Galore...");
		LOGGER.info("Log files can be found in "+outFile+".stdOut and "+outFile+".stdErr");
		//super.executeCommand(new String[]{"/storageNGS/ngs1/software/trim_galore_v0.4.3/trim_galore --path_to_cutadapt /home/ibis/paul.hager/.local/bin/cutadapt -o /home/ibis/paul.hager/KNIME4NGS/trimG_TEST -q 30 --stringency 10 -e 0.4 --gzip --length 10 --fastqc_args=\"-o /home/ibis/paul.hager/KNIME4NGS/fastqcAFTER_TEST\" --max_n 80 --no_report_file  --paired /home/ibis/paul.hager/knime-workspace/KNIME4NGS_Test_VarCalling/knime4ngs-resource/VarCalling/NA12877_R2.fastq.gz /home/ibis/paul.hager/knime-workspace/KNIME4NGS_Test_VarCalling/knime4ngs-resource/VarCalling/NA12877_R1.fastq.gz"},outFile1, exec, new String[] {path_variable}, lockFile, outFile1+".stdOut", outFile1+".stdErr", null, null, null);
		super.executeCommand(command, outFile, exec, new String[] {path_variable}, lockFile, outFile+".stdOut", outFile+".stdErr", null, null, null);
    	
		
		BufferedDataContainer cont;
	    FileCell[] c;
	    	
	   
	    cont = exec.createDataContainer(createSpecs());
	        
	    
	    if(readType.equals("single-end")){
	       	c = new FileCell[]{
	       			FileCellFactory.create(outFile)};
	    } else {
	    	c = new FileCell[]{
	       			FileCellFactory.create(outFile),
	       			FileCellFactory.create(outFile2)};
	    }
	    cont.addRowToTable(new DefaultRow("Row0", c));
	  
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	//deleteZipFiles(inFile1, inFile2);

		return new BufferedDataTable[]{outTable};
    }
    
    /**
     * Creates correct trimmed fastq output file name
     * @param file
     *  
     */
    private String createOutputName(String file){
    	
    	String outFolderTrimG = m_outfolder_trimg.getStringValue();
    	
    	if(!outFolderTrimG.equals("")){
    		file = outFolderTrimG+System.getProperty("file.separator")+file.substring(file.lastIndexOf(System.getProperty("file.separator"))+1, file.length());
    	}
    	    	
    	if(readType.equals("single-end")){
    		if(file.endsWith(".fastq.gz")){
        		file = file.replace(".fastq.gz", "_trimmed.fq.gz");
        	} else if (file.endsWith(".fastq")){
        		file = file.replace(".fastq", "_trimmed.fq");
        	} else if (file.endsWith("fq.gz")){
        		file = file.replace(".fq.gz", "_trimmed.fq.gz");
        	} else if (file.endsWith(".fq")){
        		file = file.replace(".fq", "_trimmed.fq");
        	}
    	} else {
    		if(file.endsWith(".fastq.gz")){
    			file = file.replace(".fastq.gz", "_val_"+count+".fq.gz");
        	} else if (file.endsWith(".fastq")){
        		file = file.replace(".fastq", "_val_"+count+".fq.gz");
        	} else if (file.endsWith("fq.gz")){
        		file = file.replace(".fq.gz", "_val_"+count+".fq.gz");
        	} else if (file.endsWith(".fq")){
        		file = file.replace(".fq", "_val_"+count+".fq.gz");
        	}
    		count--;
    	}
    	
    	
    	return file;
    }
    
    

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
           	
    	CompatibilityChecker CC = new CompatibilityChecker();
    	readType = CC.getReadType(inSpecs, 0);
    	if(CC.getWarningStatus()){
    		setWarningMessage(CC.getWarningMessages());
    	}
    	
    	
    	
    	if(!m_adapter.getStringValue().matches("[ACGT]+") && !m_adapter.getStringValue().equals("")){
    		throw new InvalidSettingsException("Invalid adapter! Adapter string may only contain A, C, G, or T");
    	}
    	
    	if(!m_adapter2.getStringValue().matches("[ACGT]+") && !m_adapter2.getStringValue().equals("")){
    		throw new InvalidSettingsException("Invalid adapter2! Adapter2 string may only contain A, C, G, or T");
    	}
    	
    	if(CompatibilityChecker.inputFileNotOk(m_trimg.getStringValue(), false)){
    		throw new InvalidSettingsException("Path to TrimGalore does not exist or is invalid!");
    	}
    	
    	/*
    	if(readType.equals("paired-end") && m_max_length.getIntValue() != 0){
    		throw new InvalidSettingsException("Maximum length filtering works currently only in single-end mode (which is more sensible for smallRNA-sequencing anyway...)");
    	}
    	*/
    	
    	super.updatePrefs();

        return new DataTableSpec[]{createSpecs()};

    }
    
    /**
     * Create Tablespecs
     * @return
     */
    private DataTableSpec createSpecs(){
    	DataTableSpec out;
    	
    	if(readType.equals("single-end")){
    		out = new DataTableSpec(
            		new DataColumnSpec[]{
            				new DataColumnSpecCreator(OUT_TRIMREAD1, FastQCell.TYPE).createSpec()});
    	} else {
    		out = new DataTableSpec(
            		new DataColumnSpec[]{
            				new DataColumnSpecCreator(OUT_TRIMREAD1, FastQCell.TYPE).createSpec(),
            				new DataColumnSpecCreator(OUT_TRIMREAD2, FastQCell.TYPE).createSpec()});
    	}
    	
    	
    	return out;
    }

}


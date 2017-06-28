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
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
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
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;


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
			new SettingsModelIntegerBounded(CFGKEY_THREADS_FASTQC, 4, 1, Integer.MAX_VALUE);
	
	private final SettingsModelIntegerBounded m_quality = 
			new SettingsModelIntegerBounded(CFGKEY_QUALITY, 20, 0, Integer.MAX_VALUE);
	
	private final SettingsModelString m_adapter = 
			new SettingsModelString(CFGKEY_ADAPTER, "");
	
	private final SettingsModelString m_adapter2 = 
			new SettingsModelString(CFGKEY_ADAPTER2,"");
	
	private final SettingsModelString m_preset_adapter = 
			new SettingsModelString(CFGKEY_PRESET_ADAPTER,"");
	
	//private final SettingsModelIntegerBounded m_max_length = 
	//		new SettingsModelIntegerBounded(CFGKEY_MAX_LENGTH, 0, 0, Integer.MAX_VALUE);
	
	private final SettingsModelIntegerBounded m_stringency = 
			new SettingsModelIntegerBounded(CFGKEY_STRINGENCY, 1, 1, Integer.MAX_VALUE);
	
	private final SettingsModelDoubleBounded m_error_rate = 
			new SettingsModelDoubleBounded(CFGKEY_ERROR_RATE, 0.1, 0, 1);
	
	private final SettingsModelBoolean m_gzip = 
			new SettingsModelBoolean(CFGKEY_GZIP, true);
	
	private final SettingsModelIntegerBounded m_length = 
			new SettingsModelIntegerBounded(CFGKEY_LENGTH, 20, 0, Integer.MAX_VALUE);
	
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
    
        // TODO one incoming port and one outgoing port is assumed
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
    	
    	if(m_cutadapt.equals("") || Files.notExists(Paths.get(m_cutadapt.getStringValue()))) {
    		setWarningMessage("Cutadapt PATH was not specified!");
    	} else {
    		cmd.add("--path_to_cutadapt="+IO.processFilePath(m_cutadapt.getStringValue()));
    	}
    	
    	if(m_fastqc.equals("")|| Files.notExists(Paths.get(m_fastqc.getStringValue()))) {
    		setWarningMessage("FastQC PATH was not specified!");
    	} else {
    		//remove 'fastqc' at end of path
    		int pos = m_fastqc.getStringValue().lastIndexOf(System.getProperty("file.separator"));
    		path_variable +=m_fastqc.getStringValue().substring(0,pos)+":";
    	}
    	
    	path_variable += System.getenv("PATH");
    	
    	if(!m_outfolder_trimg.getStringValue().equals("")){
    		cmd.add("-o="+IO.processFilePath(m_outfolder_trimg.getStringValue()));
    	} else {
    		String inFolder = inFile1.substring(0, inFile1.lastIndexOf(System.getProperty("file.separator"))+1);
    		cmd.add("-o="+inFolder);
    	}
    	
    	File lockFile = new File(outFile.substring(0,outFile.lastIndexOf(".")) + ".TrimG" + SuccessfulRunChecker.LOCK_ENDING);
    	
    	if(m_quality.getIntValue() != 20){
    		cmd.add("-q="+m_quality.getIntValue());
    	}
    	
    	if(!m_adapter.getStringValue().equals("")){
    		cmd.add("-a="+m_adapter.getStringValue());
    	}
    	
    	if(!m_adapter2.getStringValue().equals("")){
    		cmd.add("-a2="+m_adapter2.getStringValue());
    	}
    	
    	if(!m_preset_adapter.getStringValue().equals("") && m_adapter.getStringValue().equals("")){
    		if(m_preset_adapter.getStringValue().equals("illumina")){
    			cmd.add("--illumina");
    		} else if(m_preset_adapter.getStringValue().equals("nextera")){
    			cmd.add("--nextera");
    		} else if (m_preset_adapter.getStringValue().equals("small_rna")){
    			cmd.add("--small_rna");
    		}
    	}
    	
    	/*
    	if(m_max_length.getIntValue() != 0){
    		cmd.add("--max_length="+m_max_length.getIntValue());
    	}
    	*/
    	
    	if(m_stringency.getIntValue() != 1){
    		cmd.add("--stringency="+m_stringency.getIntValue());
    	}
    	
    	if(m_error_rate.getDoubleValue() != 0.1){
    		cmd.add("-e="+m_error_rate.getDoubleValue());
    	}
    	
    	if(m_gzip.getBooleanValue()){
    		cmd.add("--gzip");
    	} else {
    		cmd.add("--dont_gzip");
    	}
    	
    	if(m_length.getIntValue() != 20){
    		cmd.add("--length="+m_length.getIntValue());
    	}
    	
    	if(m_fastqc_enabled.getBooleanValue()){
    		String fastqcString = "";
    		if(m_threads.getIntValue() != 1 || !m_outfolder_fastqc.getStringValue().equals("")){
    			fastqcString += "--fastqc_args=";
    			if(m_threads.getIntValue() != 1){
    				fastqcString += "-t="+m_threads.getIntValue()+" ";
    			}
    			if(!m_outfolder_fastqc.getStringValue().equals("")){
    				fastqcString += "-o="+IO.processFilePath(m_outfolder_fastqc.getStringValue())+" ";
    			}
    			if(!m_fastqc_additional_options.getStringValue().equals("")){
    				setWarningMessage("NOTE! All additional FastQC options MUST use equals operator to bind arguments to options. i.e. --k=5");
    				fastqcString += m_fastqc_additional_options.getStringValue();	
    			}
    			fastqcString = fastqcString.trim();
    			fastqcString += "";
    		} else {
    			fastqcString += "--fastqc";
    		}
    		cmd.add(fastqcString);
    	}
    	
    	if(!m_additional_options.getStringValue().equals("")){
    		setWarningMessage("NOTE! All additional TrimGalore options MUST use equals operator to bind arguments to options. i.e. --length_1=30");
    		String[] addOpts = m_additional_options.getStringValue().split(" ");
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
    	if(!m_outfolder_trimg.getStringValue().equals("")){
    		file = m_outfolder_trimg.getStringValue()+System.getProperty("file.separator")+file.substring(file.lastIndexOf(System.getProperty("file.separator"))+1, file.length());
    	}
    	
    	System.out.println();
    	
    	if(readType.equals("single-end")){
    		if(file.endsWith(".fastq.gz")){
        		file = file.replace(".fastq.gz", "_trimmed.fq.gz");
        	} else {
        		file = file.replace(".fastq", "_trimmed.fq");
        	}
    	} else {
    		if(file.endsWith(".fastq.gz")){
    			file = file.replace(".fastq.gz", "_val_"+count+".fq.gz");
        	} else {
        		file = file.replace(".fastq", "_val_"+count+".fq");
        	}
    	}
    	
    	count--;
    	
    	return file;
    }
    
    

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO Code executed on reset.
        // Models build during execute are cleared here.
        // Also data handled in load/saveInternals will be erased here.
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
        
        // TODO: check if user settings are available, fit to the incoming
        // table structure, and the incoming types are feasible for the node
        // to execute. If the node can execute in its current state return
        // the spec of its output data table(s) (if you can, otherwise an array
        // with null elements), or throw an exception with a useful user message
    	
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

        // TODO: generated method stub
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


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
package de.helmholtz_muenchen.ibis.ngs.picard;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.SAMCell;
import de.helmholtz_muenchen.ibis.utils.lofs.PathProcessor;



/**
 * This is the model implementation of PicardTools.
 * 
 *
 * @author Marie-Sophie Friedl
 * @author Maximilian Hastreiter
 */
public class PicardToolsNodeModel extends HTExecutorNodeModel {
    
    // the logger instance
    protected static final NodeLogger logger = NodeLogger.getLogger(PicardToolsNodeModel.class);
    
    
    //general options
    static final String CFGKEY_PICARD="picard";
    private final SettingsModelString m_picard = new SettingsModelString(CFGKEY_PICARD, "");
    static final String CFGKEY_PICARD_MEM="picard_mem";
    private final SettingsModelIntegerBounded m_picard_mem = new SettingsModelIntegerBounded(CFGKEY_PICARD_MEM, 8, 1, Integer.MAX_VALUE);
    static final String CFGKEY_PICARD_TMP="picard_tmp";
    private final SettingsModelOptionalString m_picard_tmp = new SettingsModelOptionalString(CFGKEY_PICARD_TMP, "",false);
    static final String CFGKEY_REFGENOME="refgenome";
    private final SettingsModelString m_refgenome = new SettingsModelString(CFGKEY_REFGENOME, "");
    
    //selected tool
    static final String CFGKEY_PTOOL="ptool";
    static final String[] TOOLS_AVAILABLE={"AddOrReplaceReadGroups", "CollectInsertSizeMetrics", "MarkDuplicates", "SortSam"};
    private final SettingsModelString m_ptool = new SettingsModelString(CFGKEY_PTOOL, "");
    //output format bam/sam
    static final String CFGKEY_BSFORMAT="bsformat";
    static final String DEF_BSFORMAT="bam";
    private final SettingsModelString m_bsformat = new SettingsModelString(CFGKEY_BSFORMAT, DEF_BSFORMAT);
    //create index?
    static final String CFGKEY_INDEX="index";
    static final boolean DEF_INDEX=true;
    private final SettingsModelBoolean m_index = new SettingsModelBoolean(CFGKEY_INDEX, DEF_INDEX);
    //validation stringency
    static final String CFGKEY_VALSTRING="valstring";
    static final String[] VALSTRING_AVAILABLE={"SILENT", "LENIENT", "STRICT"};
    static final String DEF_VALSTRING=VALSTRING_AVAILABLE[0];
    private final SettingsModelString m_valstring = new SettingsModelString(CFGKEY_VALSTRING, DEF_VALSTRING);
    
    // add or replace read groups
    
    //use file name to create id, library, sample name?
    static final String CFGKEY_USE_FILE_NAME="use_file_name";
    static final boolean DEF_USE_FILE_NAME=false;
    private final SettingsModelBoolean m_use_file_name = new SettingsModelBoolean(CFGKEY_USE_FILE_NAME, DEF_USE_FILE_NAME);
    //id name
    static final String CFGKEY_ID_NAME="id_name";
    static final String DEF_ID_NAME="id";
    private final SettingsModelString m_id_name = new SettingsModelString(CFGKEY_ID_NAME, DEF_ID_NAME);
    //library name
    static final String CFGKEY_LIBRARY_NAME="library_name";
    static final String DEF_LIBRARY_NAME="library";
    private final SettingsModelString m_library_name= new SettingsModelString(CFGKEY_LIBRARY_NAME, DEF_LIBRARY_NAME);
    //sample name
    static final String CFGKEY_SAMPLE_NAME="sample_name";
    static final String DEF_SAMPLE_NAME="sample";
    final SettingsModelString m_sample_name = new SettingsModelString(CFGKEY_SAMPLE_NAME, DEF_SAMPLE_NAME);
    //sequencing platform unit
    static final String CFGKEY_PLATFROM_UNIT="platform_unit";
    static final String DEF_PLATFORM_UNIT="unit";
    private final SettingsModelString m_platform_unit = new SettingsModelString(CFGKEY_PLATFROM_UNIT, DEF_PLATFORM_UNIT);
    //sequencing platform
    static final String CFGKEY_PLATFORM="platform";
    static final String DEF_PLATFROM="ILLUMINA";
    private final SettingsModelString m_platform = new SettingsModelString(CFGKEY_PLATFORM, DEF_PLATFROM);
    
    //insert size metrics
    
    //accumulation level
    static final String CFGKEY_ACC_LEVEL="acc_level";
    static final String DEF_ACC_LEVEL="ALL_READS";
    private final SettingsModelString m_acc_level = new SettingsModelString(CFGKEY_ACC_LEVEL, DEF_ACC_LEVEL);
    //is sorted?
    static final String CFGKEY_ASS_SORTED_SM = "ass_sorted_sm";
    static final boolean DEF_ASS_SORTED_SM=false;
    private final SettingsModelBoolean m_ass_sorted_sm=new SettingsModelBoolean(CFGKEY_ASS_SORTED_SM, DEF_ASS_SORTED_SM);
    //minimum percentage of reads
    static final String CFGKEY_MIN_PCT="min_pct";
    static final double DEF_MIN_PCT=0.05;
    static final double MIN_MIN_PC=0;
    static final double MAX_MIN_PC=0.5;
    private final SettingsModelDoubleBounded m_min_pct = new SettingsModelDoubleBounded(CFGKEY_MIN_PCT, DEF_MIN_PCT, MIN_MIN_PC, MAX_MIN_PC);
    //deviation for plots and calculation
    static final String CFGKEY_DEVIATION="deviation";
    static final double DEF_DEVIATION=10.0;
    static final double MIN_DEVIATION=0;
    static final double MAX_DEVIATION=Double.MAX_VALUE;
    private final SettingsModelDoubleBounded m_deviation = new SettingsModelDoubleBounded(CFGKEY_DEVIATION, DEF_DEVIATION, MIN_DEVIATION, MAX_DEVIATION);
    
    //mark duplicates
    
    //remove duplicates?
    static final String CFGKEY_REMOVE_DUPL="remove_dupl";
    static final boolean DEF_REMOVE_DUPL=false;
    private final SettingsModelBoolean m_remove_dupl = new SettingsModelBoolean(CFGKEY_REMOVE_DUPL, DEF_REMOVE_DUPL);
    //is sorted?
    static final String CFGKEY_ASS_SORTED_RD="ass_sorted";
    static final boolean DEF_ASS_SORTED_RD=false;
    private final SettingsModelBoolean m_ass_sorted_rd = new SettingsModelBoolean(CFGKEY_ASS_SORTED_RD, DEF_ASS_SORTED_RD);
    
    //sorting options
    
    //sort order
    static final String CFGKEY_SORT_ORDER="sort_order";
    static final String[] SORT_AVAILABLE={"coordinate","queryname","unsorted"};
    static final String DEF_SORT_ORDER=SORT_AVAILABLE[0];
    private final SettingsModelString m_sort_order=new SettingsModelString(CFGKEY_SORT_ORDER, DEF_SORT_ORDER);
    
    private int posSamBam;
    private String picard_bin, ref_genome;

    protected PicardToolsNodeModel() {
    
        super(1, 1);
        
    	//general options
        addSetting(m_picard);
        addSetting(m_picard_mem);
        addSetting(m_picard_tmp);
        addSetting(m_refgenome);
        addSetting(m_ptool);
        addSetting(m_bsformat);
        addSetting(m_index);
        addSetting(m_valstring);
        //add or replace read groups
        addSetting(m_use_file_name);
        addSetting(m_id_name);
        addSetting(m_library_name);
        addSetting(m_sample_name);
        addSetting(m_platform_unit);
        addSetting(m_platform);
        //insert size metrics
        addSetting(m_acc_level);
        addSetting(m_ass_sorted_sm);
        addSetting(m_min_pct);
        addSetting(m_deviation);
        //mark duplicates
        addSetting(m_remove_dupl);
        addSetting(m_ass_sorted_rd);
        //sorting
        addSetting(m_sort_order);
        
    	addPrefPageSetting(m_picard, IBISKNIMENodesPlugin.PICARD);
    	addPrefPageSetting(m_refgenome, IBISKNIMENodesPlugin.REF_GENOME);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	
    	//Check input table integrity
    	CompatibilityChecker.inDataCheck(inData);
        
        //retrieves BAM/SAM file and reference sequence file from table of previous node
        DataRow r=inData[0].iterator().next();
        String inputfile=r.getCell(posSamBam).toString();
        
        //checks if all files are still available
        if(!Files.exists(Paths.get(inputfile))){ 
        	throw new Exception("file: "+ inputfile+" does not exist");
        }
        
        
        logger.info("Using file "+inputfile+" with reference "+ref_genome);
        
        
        //process path to input file -> location and base name of output file        
        String fileextension = PathProcessor.getExt(inputfile);
        String base = PathProcessor.getBase(inputfile);
        String basename=Paths.get(base).getFileName().toString();
        
        
        //checks file extension
        if(!fileextension.equals("sam") && !fileextension.equals("bam")){
        	throw new Exception("Input file is not in sam or bam format!");
        }
        
        String tool=m_ptool.getStringValue();
        
        //check if tool is specified and node has been properly configured
        if(tool.equals("")){
        	throw new Exception("You have to configure this node before running it!!!");
        }
        logger.info("Starting tool "+tool+"...");
        
        
        //launch collect insert size metrics
        if(tool.equals(TOOLS_AVAILABLE[1])){
        	
        	String output_hist=PathProcessor.createOutputFile(base, "pdf", "ismhist");
        	String output_data=PathProcessor.createOutputFile(base, "txt", "ismetrics");
        	
        	runMetrics(exec,inputfile, ref_genome, output_data, output_hist, m_valstring.getStringValue(), m_index.getBooleanValue(), m_deviation.getDoubleValue(), m_min_pct.getDoubleValue(), m_acc_level.getStringValue(), m_ass_sorted_sm.getBooleanValue());

		    //create output table with output sam/bam file and reference sequence AND metrics files
		    
		    //definition of the column labels and the column data types
		    DataColumnSpec[] colspec= new DataColumnSpec[2];
		    
		    if(m_bsformat.getStringValue().equals("bam")){
		    	colspec[0]=new DataColumnSpecCreator("Path2BAMFile", BAMCell.TYPE).createSpec();
		    }
		    else{
		    	colspec[0]=new DataColumnSpecCreator("Path2SAMFile", SAMCell.TYPE).createSpec();
		    }
		    colspec[1]=new DataColumnSpecCreator("Path2ISMetrics", FileCell.TYPE).createSpec();
		    DataTableSpec outspec=new DataTableSpec(colspec);
		    
		    //container which is filled with cells and rows
		    BufferedDataContainer c = exec.createDataContainer(outspec);
		    		    
		    FileCell sc1 = FileCellFactory.create(inputfile);
		    FileCell sc3 = FileCellFactory.create(output_data);
		    
		    //create row and add it to the container
		    DefaultRow row = new DefaultRow("row0", new FileCell[]{sc1, sc3});
		    c.addRowToTable(row);
		    
		    //create final table
		    c.close();
		    BufferedDataTable out=c.getTable();
		    
		    return new BufferedDataTable[]{out};
        }
        
        else{
            String output="";
        	
		    //launch add or replace read groups
		    if(tool.equals(TOOLS_AVAILABLE[0])){
		    	
		    	output=PathProcessor.createOutputFile(base, m_bsformat.getStringValue(), "rg");
		    	
		    	if(m_use_file_name.getBooleanValue()){
		    		runRG(exec,inputfile, output, m_valstring.getStringValue(), m_index.getBooleanValue(), basename, basename, basename, basename, m_platform.getStringValue());
		    	}
		    	else{
		    		runRG(exec,inputfile, output, m_valstring.getStringValue(), m_index.getBooleanValue(), m_id_name.getStringValue(), m_library_name.getStringValue(), m_sample_name.getStringValue(), m_platform_unit.getStringValue(), m_platform.getStringValue());		    		
		    	}
		    }
		    

		    //launch mark duplicates
		    else if(tool.equals(TOOLS_AVAILABLE[2])){
		    	
		    	output=PathProcessor.createOutputFile(base, m_bsformat.getStringValue(), "marked");
		    	String output_metrics=PathProcessor.createOutputFile(base, "txt", "marked.metrics");
		    	
		    	runDupl(exec,inputfile, output, output_metrics, m_valstring.getStringValue(), m_index.getBooleanValue(), m_remove_dupl.getBooleanValue(), m_ass_sorted_rd.getBooleanValue());
		    }
		    
		    //launch sort sam
		    else if(tool.equals(TOOLS_AVAILABLE[3])){
		    	
		    	output=PathProcessor.createOutputFile(base, m_bsformat.getStringValue(), "sorted");
		    	runSortSam(exec,inputfile, output, m_valstring.getStringValue(), m_index.getBooleanValue(), m_sort_order.getStringValue());
		    }
		    
		    //create output table with output sam/bam file and reference sequence
		    
		    //definition of the column labels and the column data types
		    DataColumnSpec[] colspec= new DataColumnSpec[1];
		    
		    if(m_bsformat.getStringValue().equals("bam")){
		    	colspec[0]=new DataColumnSpecCreator("Path2BAMFile", BAMCell.TYPE).createSpec();
		    }
		    else{
		    	colspec[0]=new DataColumnSpecCreator("Path2SAMFile", SAMCell.TYPE).createSpec();
		    }
		    
		    //create column specifications for table
		    DataTableSpec outspec=new DataTableSpec(colspec);
		    
		    //container which is filled with cells and rows
		    BufferedDataContainer c = exec.createDataContainer(outspec);
		    
		    //single cells containing paths to use for the next node
		    FileCell sc1= FileCellFactory.create(output);
		    
		    //create row and add it to the container
		    DefaultRow row = new DefaultRow("row0", new FileCell[]{sc1});
		    c.addRowToTable(row);
		    
		    //create final table
		    c.close();
		    BufferedDataTable out=c.getTable();
		    
		    return new BufferedDataTable[]{out};
		    
        }

    }

	/**
	 * {@inheritDoc}
	 */
	// executed when in port is connected to executed node or when previous node
	// (already connected) is executed
	@Override
	// executed when in port is connected
	protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {

		super.updatePrefs();
		picard_bin = IO.processFilePath(m_picard.getStringValue());
		ref_genome = IO.processFilePath(m_refgenome.getStringValue());
		
		if (CompatibilityChecker.inputFileNotOk(picard_bin, false)) {
			throw new InvalidSettingsException("Set path to picard.jar!");
		}
		
        if(!Files.exists(Paths.get(ref_genome))){
        	throw new InvalidSettingsException("file: "+ ref_genome+" does not exist");
        }

		// names of the input table columns
		String[] cols = inSpecs[0].getColumnNames();

		// checking for input bam/sam file
		if (!CompatibilityChecker.checkInputCellType(inSpecs[0], "SAMCell")
				&& !CompatibilityChecker.checkInputCellType(inSpecs[0], "BAMCell")) {
			throw new InvalidSettingsException("Previous node is incompatible! Missing path to sam/bam file!");
		}

		// determining position of reference and sam/bam file
		for (int i = 0; i < cols.length; i++) {
			if (cols[i].equals("Path2BAMFile") || cols[i].equals("Path2SAMFile")) {
				posSamBam = i;
			}
		}

		boolean add_col = m_ptool.getStringValue().equals(TOOLS_AVAILABLE[1]);

		DataColumnSpec[] colspec;
		if (add_col) {
			colspec = new DataColumnSpec[2];
		} else {
			colspec = new DataColumnSpec[1];
		}

		if (m_bsformat.getStringValue().equals("bam")) {
			colspec[0] = new DataColumnSpecCreator("Path2BAMFile", BAMCell.TYPE).createSpec();
		} else {
			colspec[0] = new DataColumnSpecCreator("Path2SAMFile", SAMCell.TYPE).createSpec();
		}

		// create column specifications for table
		if (add_col) {
			colspec[1] = new DataColumnSpecCreator("Path2ISMetrics", FileCell.TYPE).createSpec();
		}

		return new DataTableSpec[] { new DataTableSpec(colspec) };
	}

        
	//run CollectInsertSizeMetrics
	protected void runMetrics (ExecutionContext exec,String input, String ref, String outputm, String outputh, String val, boolean index, double dev, double pct, String acc, boolean sorted ) throws Exception{
		
		/* command line options
		 * INPUT
		 * CREATE_INDEX
		 * VALIDATION_STRINGENCY
		 * HISTOGRAM_FILE
		 * DEVIATIONS
		 * MINIMUM_PCT
		 * METRIC_ACCUMULATION_LEVEL
		 * ASSUME_SORTED
		 * OUTPUT
		 * REFERENCE_SEQUENCE
		 */

		String method = "CollectInsertSizeMetrics";
		String[] args = new String[11];
		args[0]="INPUT="+input;
		args[1]="OUTPUT="+outputm;
		args[2]="HISTOGRAM_FILE="+outputh;
		args[3]="REFERENCE_SEQUENCE="+ref;
		args[4]="VALIDATION_STRINGENCY="+val;
		args[5]="CREATE_INDEX="+index;
		args[6]="DEVIATIONS="+dev;
		args[7]="MINIMUM_PCT="+pct;
		args[8]="METRIC_ACCUMULATION_LEVEL="+acc;
		args[9]="ASSUME_SORTED="+sorted;
		
		if(m_picard_tmp.isActive()){
			args[10]="TMP_DIR="+m_picard_tmp.getStringValue();
		}else{
			args[10]="";
		}
		
		
		
		File lockFile = new File(outputm+SuccessfulRunChecker.LOCK_ENDING);
		String command = "java -Xmx"+m_picard_mem.getIntValue()+"G -jar "+picard_bin+" "+method+" "+args[0]+" "+args[1]+" "+args[2]+" "+args[3]+" "+args[4]+" "+args[5]+" "+args[6]+" "+args[7]+" "+args[8]+" "+args[9]+" "+args[10];
		super.executeCommand(new String[] { command }, outputm, exec, lockFile, outputm+".stdOut",outputm+".stdErr");
		
		PicardToolsNodeModel.logger.info("CollectInsertSizeMetrics finished successfully");
	}
	
	//run AddOrReplaceReadGroups
	protected void runRG(ExecutionContext exec, String input, String output, String val, boolean index, String id, String library, String sample, String unit, String platform) throws Exception{
		
		/* command line arguments
		 * INPUT
		 * OUTPUT
		 * VALIDATION_STRINGENCY
		 * CREATE_INDEX
		 * RGID id
		 * RGLB library
		 * RGPL platform
		 * RGPU platform unit
		 * RGSM sample
		 */
		
		if(id.equals("")) {
			id = "id";
			PicardToolsNodeModel.logger.warn("ID has to be specified (now set to 'id')");
		}
		if(library.equals("")) {
			library = "library";
			PicardToolsNodeModel.logger.warn("Library has to be specified (now set to 'library')");
		}
		if(sample.equals("")) {
			sample="sample";
			PicardToolsNodeModel.logger.warn("Sample has to be specified (now set to 'sample')");
		}
		if(unit.equals("")) {
			unit="unit";
			PicardToolsNodeModel.logger.warn("Platform unit has to be specified (now set to 'unit')");
		}
		if(platform.equals("")) {
			platform = "platform";
			PicardToolsNodeModel.logger.warn("Platform has to be specified (now set to 'platform')");
		}
		
		String method = "AddOrReplaceReadGroups";
		
		String [] args = new String[10];
		args[0]="INPUT="+input;
		args[1]="OUTPUT="+output;
		args[2]="VALIDATION_STRINGENCY="+val;
		args[3]="CREATE_INDEX="+index;
		args[4]="RGID="+id;
		args[5]="RGLB="+library;
		args[6]="RGSM="+sample;
		args[7]="RGPU="+unit;
		args[8]="RGPL="+platform;
		
		if(m_picard_tmp.isActive()){
			args[9]="TMP_DIR="+m_picard_tmp.getStringValue();
		}else{
			args[9]="";
		}
		
		File lockFile = new File(output+SuccessfulRunChecker.LOCK_ENDING);
		String command = "java -Xmx"+m_picard_mem.getIntValue()+"G -jar "+picard_bin+" "+method+" "+args[0]+" "+args[1]+" "+args[2]+" "+args[3]+" "+args[4]+" "+args[5]+" "+args[6]+" "+args[7]+" "+args[8]+" "+args[9];

		super.executeCommand(new String[] { command }, output, exec, lockFile, output+".stdOut",output+".stdErr");

		PicardToolsNodeModel.logger.info("AddOrReplaceReadGroups finished successfully");
	}

	
	//run MarkDuplicates
	protected void runDupl(ExecutionContext exec, String input, String output, String metrics, String val, boolean index, boolean rmdupl, boolean ass_sort) throws Exception {
		
		/* command line options for MarkDuplicates
		 * INPUT
		 * OUTPUT
		 * METRICS_FILE
		 * REMOVE_DUPLICATES
		 * ASSUME_SORTED
		 * VALIDATION_STRINGENCY
		 * CREATE_INDEX
		 */
		String method = "MarkDuplicates";
		String [] args = new String[8];
		args[0]="INPUT="+input;
		args[1]="OUTPUT="+output;
		args[2]="METRICS_FILE="+metrics;
		args[3]="VALIDATION_STRINGENCY="+val;
		args[4]="CREATE_INDEX="+index;
		args[5]="REMOVE_DUPLICATES="+rmdupl;
		args[6]="ASSUME_SORTED="+ass_sort;
		if(m_picard_tmp.isActive()){
			args[7]="TMP_DIR="+m_picard_tmp.getStringValue();
		}else{
			args[7]="";
		}
		
		File lockFile = new File(output+SuccessfulRunChecker.LOCK_ENDING);
		String command = "java -Xmx"+m_picard_mem.getIntValue()+"G -jar "+picard_bin+" "+method+" "+args[0]+" "+args[1]+" "+args[2]+" "+args[3]+" "+args[4]+" "+args[5]+" "+args[6]+" "+args[7];

		super.executeCommand(new String[] { command }, output, exec, lockFile, output+".stdOut",output+".stdErr");
		
		PicardToolsNodeModel.logger.info("MarkDupplicates finished successfully");
	}
	
	//run SortSam
	protected void runSortSam(ExecutionContext exec, String input, String output, String val, boolean index, String order) throws Exception {
		
		/* command line options for SortSam
		 *INPUT
		 *OUTPUT
		 *SORT_ORDER
		 *VALIDATION_STRINGENCY
		 *CREATE_INDEX
		 */
		String method = "SortSam";
		
		String[] args= new String[6];
		args[0]="INPUT="+input;
		args[1]="OUTPUT="+output;
		args[2]="VALIDATION_STRINGENCY="+val;
		args[3]="CREATE_INDEX="+index;
		args[4]="SORT_ORDER="+order;
		if(m_picard_tmp.isActive()){
			args[5]="TMP_DIR="+m_picard_tmp.getStringValue();
		}else{
			args[5]="";
		}
		
		File lockFile = new File(output+SuccessfulRunChecker.LOCK_ENDING);
		String command = "java -Xmx"+m_picard_mem.getIntValue()+"G -jar "+picard_bin+" "+method+" "+args[0]+" "+args[1]+" "+args[2]+" "+args[3]+" "+args[4]+" "+args[5];
		
		super.executeCommand(new String[] { command }, output, exec, lockFile, output+".stdOut",output+".stdErr");

		PicardToolsNodeModel.logger.info("SortSam finished successfully");
	}   
}


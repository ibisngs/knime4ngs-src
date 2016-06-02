package de.helmholtz_muenchen.ibis.ngs.picard;

import java.io.File;
import java.io.IOException;
import java.nio.file.*;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
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
 * @author 
 */
public class PicardToolsNodeModel extends HTExecutorNodeModel {
    
    // the logger instance
    protected static final NodeLogger logger = NodeLogger.getLogger(PicardToolsNodeModel.class);
    
    
    //general options
    static final String CFGKEY_PICARD="picard";
    private final SettingsModelString m_picard = new SettingsModelString(CFGKEY_PICARD, "");
    static final String CFGKEY_PICARD_MEM="picard_mem";
    private final SettingsModelIntegerBounded m_picard_mem = new SettingsModelIntegerBounded(CFGKEY_PICARD_MEM, 8, 1, Integer.MAX_VALUE);
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

    protected PicardToolsNodeModel() {
    
        super(1, 1);
        
    	//general options
        addSetting(m_picard);
        addSetting(m_picard_mem);
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
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

        
        //retrieves BAM/SAM file and reference sequence file from table of previous node
        DataRow r=inData[0].iterator().next();
        String inputfile=r.getCell(posSamBam).toString();
        String reffile=m_refgenome.getStringValue();
        
        
        //checks if all files are still available
        if(!Files.exists(Paths.get(inputfile))){ 
        	throw new Exception("file: "+ inputfile+" does not exist");
        }
        
        if(!Files.exists(Paths.get(reffile))){
        	throw new Exception("file: "+ reffile+" does not exist");
        }
        
        logger.info("Using file "+inputfile+" with reference "+reffile);
        
        
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
        	
        	runMetrics(exec,inputfile, reffile, output_data, output_hist, m_valstring.getStringValue(), m_index.getBooleanValue(), m_deviation.getDoubleValue(), m_min_pct.getDoubleValue(), m_acc_level.getStringValue(), m_ass_sorted_sm.getBooleanValue());

		    //create output table with output sam/bam file and reference sequence AND metrics files
		    
		    //definition of the column labels and the column data types
		    DataColumnSpec[] colspec= new DataColumnSpec[3];
		    
		    if(m_bsformat.getStringValue().equals("bam")){
		    	colspec[0]=new DataColumnSpecCreator("Path2BAMFile", BAMCell.TYPE).createSpec();
		    }
		    else{
		    	colspec[0]=new DataColumnSpecCreator("Path2SAMFile", SAMCell.TYPE).createSpec();
		    }
		    colspec[1]=new DataColumnSpecCreator("Path2SEQFile", FileCell.TYPE).createSpec();
		    colspec[2]=new DataColumnSpecCreator("Path2ISMetrics", FileCell.TYPE).createSpec();
		    DataTableSpec outspec=new DataTableSpec(colspec);
		    
		    //container which is filled with cells and rows
		    BufferedDataContainer c = exec.createDataContainer(outspec);
		    		    
		    FileCell sc1 = (FileCell)FileCellFactory.create(inputfile);
		    FileCell sc2 = (FileCell)FileCellFactory.create(reffile);
		    FileCell sc3 = (FileCell)FileCellFactory.create(output_data);
		    
		    //create row and add it to the container
		    DefaultRow row = new DefaultRow("row0", new FileCell[]{sc1, sc2, sc3});
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
		    DataColumnSpec[] colspec= new DataColumnSpec[2];
		    
		    if(m_bsformat.getStringValue().equals("bam")){
		    	colspec[0]=new DataColumnSpecCreator("Path2BAMFile", BAMCell.TYPE).createSpec();
		    }
		    else{
		    	colspec[0]=new DataColumnSpecCreator("Path2SAMFile", SAMCell.TYPE).createSpec();
		    }
		    
		    //create column specifications for table
		    colspec[1]=new DataColumnSpecCreator("Path2SEQFile", FileCell.TYPE).createSpec();
		    DataTableSpec outspec=new DataTableSpec(colspec);
		    
		    //container which is filled with cells and rows
		    BufferedDataContainer c = exec.createDataContainer(outspec);
		    
		    //single cells containing paths to use for the next node
		    FileCell sc1= (FileCell)FileCellFactory.create(output);
		    FileCell sc2= (FileCell)FileCellFactory.create(reffile);
		    
		    //create row and add it to the container
		    DefaultRow row = new DefaultRow("row0", new FileCell[]{sc1, sc2});
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
    @Override
    protected void reset() {

    }

    /**
     * {@inheritDoc}
     */
    //executed when in port is connected to executed node or when previous node (already connected) is executed
    @Override
    //executed when in port is connected
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	    		
    		//names of the input table columns
			String [] cols=inSpecs[0].getColumnNames();
			
			// checking for input bam/sam file
			if(!CompatibilityChecker.checkInputCellType(inSpecs[0],"SAMCell") && !CompatibilityChecker.checkInputCellType(inSpecs[0],"BAMCell")){
				throw new InvalidSettingsException("Previous node is incompatible! Missing path to sam/bam file!");
			}
			
			//determining position of reference and sam/bam file
			for(int i=0; i<cols.length; i++){
				if(cols[i].equals("Path2BAMFile") || cols[i].equals("Path2SAMFile")){
					posSamBam=i;
				}
			}
			
			boolean add_col = m_ptool.getStringValue().equals(TOOLS_AVAILABLE[1]);
			
			
			DataColumnSpec[] colspec;
			if(add_col) {
				colspec = new DataColumnSpec[3];
			} else {
				colspec = new DataColumnSpec[2];
			} 
			
		    if(m_bsformat.getStringValue().equals("bam")){
		    	colspec[0]=new DataColumnSpecCreator("Path2BAMFile", BAMCell.TYPE).createSpec();
		    }
		    else{
		    	colspec[0]=new DataColumnSpecCreator("Path2SAMFile", SAMCell.TYPE).createSpec();
		    }
		    
		    //create column specifications for table
		    colspec[1]=new DataColumnSpecCreator("Path2SEQFile", FileCell.TYPE).createSpec();
		    if(add_col) {
		    	colspec[2]=new DataColumnSpecCreator("Path2ISMetrics", FileCell.TYPE).createSpec();
		    }

        return new DataTableSpec[]{new DataTableSpec(colspec)};
    }

    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {


    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {

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

		String picard = m_picard.getStringValue();
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
		args[10]="TMP_DIR="+Paths.get(input).getParent().toString();
		
		File lockFile = new File(outputm+SuccessfulRunChecker.LOCK_ENDING);
		String command = "java -Xmx"+m_picard_mem.getIntValue()+"G -jar "+picard+" "+method+" "+args[0]+" "+args[1]+" "+args[2]+" "+args[3]+" "+args[4]+" "+args[5]+" "+args[6]+" "+args[7]+" "+args[8]+" "+args[9];
		super.executeCommand(new String[] { command }, exec, lockFile, outputm+".stdOut",outputm+".stdErr");
		
		
//		PicardToolsNodeModel.logger.info(lockCommand);
		
//		boolean b = SuccessfulRunChecker.hasTerminatedSuccessfully(lockFile, lockCommand);
		
//		if(b) {
//			PicardToolsNodeModel.logger.info("According to klock CollectInsertSizeMetrics has been finished successfully!");
//			return;
//		}
			
//		SuccessfulRunChecker checker = new SuccessfulRunChecker(lockFile, lockCommand);
		
		//redirect error stream
//		redirecterr(outputm);
//		System.err.println("Output of CollectInsertSizeMetrics:");
		
//		Exception exception=null;
//		
//		int exitcode =0;
//		
//		try{
			//run tool
//			exitcode=new CollectInsertSizeMetrics().instanceMain(args);
//		}
//		catch(Exception e){
//			exception =e ;
//		}
//		//reset error stream and variable for index writing
//		finally{
//			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(false);
//			reseterr();
//		}
//		
//		//checks if exception has been caught
//		if(exception!=null){
//			throw exception;
//		}
//		
//		if(exitcode==1){
//			throw new Exception("Something went wrong while executing CollectInsertSizeMetrics, please check out the log file");
//		}
//		
//		checker.writeOK();
//		checker.finalize();
//		PicardToolsNodeModel.logger.info("CollectInsertSizeMetrics finished successfully");
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
		
		String picard = m_picard.getStringValue();
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
		args[9]="TMP_DIR="+Paths.get(input).getParent().toString();
		
		File lockFile = new File(output+SuccessfulRunChecker.LOCK_ENDING);
		String command = "java -Xmx"+m_picard_mem.getIntValue()+"G -jar "+picard+" "+method+" "+args[0]+" "+args[1]+" "+args[2]+" "+args[3]+" "+args[4]+" "+args[5]+" "+args[6]+" "+args[7]+" "+args[8];

		super.executeCommand(new String[] { command }, exec, lockFile, output+".stdOut",output+".stdErr");

		
		//		PicardToolsNodeModel.logger.info(lockCommand);
//		
//		boolean b = SuccessfulRunChecker.hasTerminatedSuccessfully(lockFile, lockCommand);
//			
//		if(b) {
//			PicardToolsNodeModel.logger.info("According to klock AddOrReplaceReadGroups has been finished successfully!");
//			return;
//		}
//			
//		SuccessfulRunChecker checker = new SuccessfulRunChecker(lockFile, lockCommand);
//			
//		//redirect error stream
//		redirecterr(output);
//		System.err.println("Output of AddOrReplaceReadGroups:");
//		
//		Exception exception=null;
//		int exitcode=0;
//		
//		try{
//			//run tool
//			exitcode=new AddOrReplaceReadGroups().instanceMain(args);
//		}
//		catch(Exception e){
//			exception =e;
//		}
//		//reset error stream, resert writing index variable
//		finally{
//			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(false);
//			reseterr();
//		}
//		
//		//checks if exception has been caught
//		if(exception!=null){
//			throw exception;
//		}
//		
//		if(exitcode==1){
//			throw new Exception("Something went wrong while executing MarkDuplicates, please check out the log file");
//		}
//		
//		checker.writeOK();
//		checker.finalize();
		
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
		
		String picard = m_picard.getStringValue();
		String method = "MarkDuplicates";
		
		String [] args = new String[8];
		args[0]="INPUT="+input;
		args[1]="OUTPUT="+output;
		args[2]="METRICS_FILE="+metrics;
		args[3]="VALIDATION_STRINGENCY="+val;
		args[4]="CREATE_INDEX="+index;
		args[5]="REMOVE_DUPLICATES="+rmdupl;
		args[6]="ASSUME_SORTED="+ass_sort;
		args[7]="TMP_DIR="+Paths.get(input).getParent().toString();
		
		File lockFile = new File(output+SuccessfulRunChecker.LOCK_ENDING);
		String command = "java -Xmx"+m_picard_mem.getIntValue()+"G -jar "+picard+" "+method+" "+args[0]+" "+args[1]+" "+args[2]+" "+args[3]+" "+args[4]+" "+args[5]+" "+args[6];

		super.executeCommand(new String[] { command }, exec, lockFile, output+".stdOut",output+".stdErr");

		
		//		PicardToolsNodeModel.logger.info(lockCommand);
//		
//		boolean b = SuccessfulRunChecker.hasTerminatedSuccessfully(lockFile, lockCommand);
//		
//		if(b) {
//			PicardToolsNodeModel.logger.info("According to klock MarkDuplicates has been finished successfully!");
//			return;
//		}
//			
//		SuccessfulRunChecker checker = new SuccessfulRunChecker(lockFile, lockCommand);
//		
//		//redirect error stream
//		redirecterr(output);
//		System.err.println("Output of MarkDuplicates:");
//		
//		Exception exception=null;
//		int exitcode=0;
//
//		try{
//			//run tool
//			exitcode = new MarkDuplicates().instanceMain(args);
//		}
//		catch(Exception e){
//			exception = e;
//		}
//		//reset error stream, reset writing index variable
//		finally{
//			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(false);
//			reseterr();
//		}
//		
//		//checks if exception has been caught
//		if(exception!=null){
//			throw exception;
//		}
//		
//		if(exitcode==1){
//			throw new Exception("Something went wrong while executing MarkDuplicates, please check out the log file");
//		}
		
		PicardToolsNodeModel.logger.info("MarkDupplicates finished successfully");
		
//		checker.writeOK();
//		checker.finalize();
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

		String picard = m_picard.getStringValue();
		String method = "SortSam";
		
		String[] args= new String[6];
		args[0]="INPUT="+input;
		args[1]="OUTPUT="+output;
		args[2]="VALIDATION_STRINGENCY="+val;
		args[3]="CREATE_INDEX="+index;
		args[4]="SORT_ORDER="+order;
		args[5]="TMP_DIR="+Paths.get(input).getParent().toString();
		
		File lockFile = new File(output+SuccessfulRunChecker.LOCK_ENDING);
		String command = "java -Xmx"+m_picard_mem.getIntValue()+"G -jar "+picard+" "+method+" "+args[0]+" "+args[1]+" "+args[2]+" "+args[3]+" "+args[4];
		
		super.executeCommand(new String[] { command }, exec, lockFile, output+".stdOut",output+".stdErr");
//		PicardToolsNodeModel.logger.info(lockCommand);
//		
//		boolean b = SuccessfulRunChecker.hasTerminatedSuccessfully(lockFile, lockCommand);
//		
//		if(b) {
//			PicardToolsNodeModel.logger.info("According to klock SortSam has been finished successfully!");
//			return;
//		}
//			
//		SuccessfulRunChecker checker = new SuccessfulRunChecker(lockFile, lockCommand);
//		
//		//redirect error stream
//		redirecterr(output);
//		System.err.println("Output of Sortsam");
//		
//		Exception exception=null;
//		int exitcode=0;
//		
//		try{
//			//run tool
//			exitcode=new SortSam().instanceMain(args);
//		}
//		catch(Exception e){
//			exception = e;
//		}
//		//reset error stream, reset index writing variable
//		finally{
//			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(false);
//			reseterr();
//		}
//		
//		//checks if exception has been caught
//		if(exception!=null){
//			throw exception;
//		}
//		if(exitcode==1){
//			throw new Exception("Something went wrong while executing SortSam, please check out the log file");
//		}
//		
		PicardToolsNodeModel.logger.info("SortSam finished successfully");
//		checker.writeOK();
//		checker.finalize();
	}
	
	
	//create file that contains output of picard tools
	//code from http://www.avajava.com/tutorials/lessons/how-do-i-redirect-standard-error-to-a-file.html
	
//	private PrintStream err=null;
//	
//	//method to redirect error stream
//	private void redirecterr(String output) throws Exception{
//		
//		//makes sure err is never overridden
//		if(err==null){
//			err = System.err;
//		}
//		PrintStream errtofile = new PrintStream( new FileOutputStream(new File(output+".log")));
//		System.setErr(errtofile);
//
//		PicardToolsNodeModel.logger.info("log file can be found in "+output+".log");
//	}
//
//	//method to reset error stream
//	private void reseterr(){
//		System.setErr(err);
//	}
    
}


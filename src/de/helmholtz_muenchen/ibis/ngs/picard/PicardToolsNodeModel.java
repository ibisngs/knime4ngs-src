package de.helmholtz_muenchen.ibis.ngs.picard;

import java.io.File;
import java.io.IOException;
import java.nio.file.*;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.SAMCell;
import de.helmholtz_muenchen.ibis.utils.lofs.PathProcessor;

//import net.sf.picard.sam.SortSam;


/**
 * This is the model implementation of PicardTools.
 * 
 *
 * @author 
 */
public class PicardToolsNodeModel extends NodeModel {
    
    // the logger instance
    protected static final NodeLogger logger = NodeLogger.getLogger(PicardToolsNodeModel.class);
    
    
    //general options
    
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
    private int posRef;

    protected PicardToolsNodeModel() {
    
        // #input ports, #output ports
        super(1, 1);
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
        String reffile=r.getCell(posRef).toString();
        
        
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
        	
        	RunPicard.runMetrics(inputfile, reffile, output_data, output_hist, m_valstring.getStringValue(), m_index.getBooleanValue(), m_deviation.getDoubleValue(), m_min_pct.getDoubleValue(), m_acc_level.getStringValue(), m_ass_sorted_sm.getBooleanValue());

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
		    		RunPicard.runRG(inputfile, output, m_valstring.getStringValue(), m_index.getBooleanValue(), basename, basename, basename, basename, m_platform.getStringValue());
		    	}
		    	else{
		    		RunPicard.runRG(inputfile, output, m_valstring.getStringValue(), m_index.getBooleanValue(), m_id_name.getStringValue(), m_library_name.getStringValue(), m_sample_name.getStringValue(), m_platform_unit.getStringValue(), m_platform.getStringValue());		    		
		    	}
		    }
		    

		    //launch mark duplicates
		    else if(tool.equals(TOOLS_AVAILABLE[2])){
		    	
		    	output=PathProcessor.createOutputFile(base, m_bsformat.getStringValue(), "marked");
		    	String output_metrics=PathProcessor.createOutputFile(base, "txt", "marked.metrics");
		    	
		    	RunPicard.runDupl(inputfile, output, output_metrics, m_valstring.getStringValue(), m_index.getBooleanValue(), m_remove_dupl.getBooleanValue(), m_ass_sorted_rd.getBooleanValue());
		    }
		    
		    //launch sort sam
		    else if(tool.equals(TOOLS_AVAILABLE[3])){
		    	
		    	output=PathProcessor.createOutputFile(base, m_bsformat.getStringValue(), "sorted");
		    	RunPicard.runSortSam(inputfile, output, m_valstring.getStringValue(), m_index.getBooleanValue(), m_sort_order.getStringValue());
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
    	
    	
    	//input port: BWA (all mappers?), SAMLoader, BAMLoader
    	//reference file: "Sequence file", "Path2SEQFile", "Path2RefFile"
    	//input sam/bam file: "Path2SAMFile", "Path2BAMFile"
    		
    		//names of the input table columns
			String [] cols=inSpecs[0].getColumnNames();
			
			// checking for input bam/sam file
			if(!CompatibilityChecker.checkInputCellType(inSpecs[0],"SAMCell") && !CompatibilityChecker.checkInputCellType(inSpecs[0],"BAMCell")){
				throw new InvalidSettingsException("Previous node is incompatible! Missing path to sam/bam file!");
			}
//			
//			//checking for reference sequence
//			if(!inSpecs[0].containsName("Sequence file") && !inSpecs[0].containsName("Path2SEQFile") && !inSpecs[0].containsName("Path2RefFile")){
//				throw new InvalidSettingsException("Previous node is incompatible! Missing path to reference sequence!");
//			}
			
			//determining position of reference and sam/bam file
			for(int i=0; i<cols.length; i++){
				if(cols[i].equals("Path2BAMFile") || cols[i].equals("Path2SAMFile")){
					posSamBam=i;
				}
				if(cols[i].equals("Sequence file") || cols[i].equals("Path2SEQFile") || cols[i].equals("Path2RefFile")){
					posRef=i;
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
    //executed when dialog opened
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
        
        //general options
        m_ptool.saveSettingsTo(settings);
        m_bsformat.saveSettingsTo(settings);
        m_index.saveSettingsTo(settings);
        m_valstring.saveSettingsTo(settings);
        
        //add or replace read group
        m_use_file_name.saveSettingsTo(settings);
        m_id_name.saveSettingsTo(settings);
        m_library_name.saveSettingsTo(settings);
        m_sample_name.saveSettingsTo(settings);
        m_platform_unit.saveSettingsTo(settings);
        m_platform.saveSettingsTo(settings);
        
        //insert size metrics
        m_acc_level.saveSettingsTo(settings);
        m_ass_sorted_sm.saveSettingsTo(settings);
        m_min_pct.saveSettingsTo(settings);
        m_deviation.saveSettingsTo(settings);
        
        //mark duplicates
        m_remove_dupl.saveSettingsTo(settings);
        m_ass_sorted_rd.saveSettingsTo(settings);
        
        //sorting
        m_sort_order.saveSettingsTo(settings);
        
    }

    /**
     * {@inheritDoc}
     */
    //executed when dialog closed
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
         
    	
    	//general options
    	m_ptool.loadSettingsFrom(settings);
    	m_bsformat.loadSettingsFrom(settings);
    	m_index.loadSettingsFrom(settings);
    	m_valstring.loadSettingsFrom(settings);
    	
    	//add or replace read group
    	m_use_file_name.loadSettingsFrom(settings);
    	m_id_name.loadSettingsFrom(settings);
    	m_library_name.loadSettingsFrom(settings);
    	m_sample_name.loadSettingsFrom(settings);
    	m_platform_unit.loadSettingsFrom(settings);
    	m_platform.loadSettingsFrom(settings);
    	
    	//insert size metrics
    	m_acc_level.loadSettingsFrom(settings);
    	m_ass_sorted_sm.loadSettingsFrom(settings);
    	m_min_pct.loadSettingsFrom(settings);
    	m_deviation.loadSettingsFrom(settings);
    	
    	//mark duplicates
    	m_remove_dupl.loadSettingsFrom(settings);
    	m_ass_sorted_rd.loadSettingsFrom(settings);
    	
    	//sorting
    	m_sort_order.loadSettingsFrom(settings);
    	

    }

    /**
     * {@inheritDoc}
     */
    //executed when dialog closed
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	
    	//general options
        m_ptool.validateSettings(settings);
        m_bsformat.validateSettings(settings);
        m_index.validateSettings(settings);
        m_valstring.validateSettings(settings);
        
        //add or replace read groups
        m_use_file_name.validateSettings(settings);
        m_id_name.validateSettings(settings);
        m_library_name.validateSettings(settings);
        m_sample_name.validateSettings(settings);
        m_platform_unit.validateSettings(settings);
        m_platform.validateSettings(settings);
        
        //insert size metrics
        m_acc_level.validateSettings(settings);
        m_ass_sorted_sm.validateSettings(settings);
        m_min_pct.validateSettings(settings);
        m_deviation.validateSettings(settings);
        
        //mark duplicates
        m_remove_dupl.validateSettings(settings);
        m_ass_sorted_rd.validateSettings(settings);
        
        //sorting
        m_sort_order.validateSettings(settings);

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

}


package de.helmholtz_muenchen.ibis.ngs.samtoolshybrid;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.RowKey;
import org.knime.core.data.container.CloseableRowIterator;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.IntCell;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.ngs.thundercall.ThunderCallNodeModel;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.lofs.PathProcessor;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;


/**
 * This is the model implementation of SamtoolsHybrid.
 * 
 *
 * @author Tanzeem Haque
 */
public class SamtoolsHybridNodeModel extends NodeModel {
    
	// the logger instance
    private static final NodeLogger logger = NodeLogger.getLogger(ThunderCallNodeModel.class);
        
    /**
	 * Config Keys
	 */
	public static final String CFGKEY_SAMTOOLS_HYBRID_PATH = "Samtools-Hybrid path";
	public static final String CFGKEY_REF_GENOME = "Reference";
	private String OUTPUT_PATH;

	/**
	 * Settings model
	 */
	
    private final SettingsModelString m_SAMTOOLS_HYBRID = new SettingsModelString(CFGKEY_SAMTOOLS_HYBRID_PATH, "");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(CFGKEY_REF_GENOME, "");

  //The Output Col Names
  	public static final String OUT_COL1_TABLE1 = "OUTFILE";
    /**
     * Constructor for the node model.
     */
    protected SamtoolsHybridNodeModel() {
    
        // TODO one incoming port and one outgoing port is assumed
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

        String samHybrid = m_SAMTOOLS_HYBRID.getStringValue();
        String refFile = m_REF_GENOME.getStringValue();
        String outFile = "";

    	String bamFile = "";

        CloseableRowIterator it = inData[0].iterator();
		while (it.hasNext()) {
			DataRow row = it.next();
			bamFile = row.getCell(0).toString();
		}
		
		this.OUTPUT_PATH = bamFile.substring(0, bamFile.lastIndexOf("/"))+"/";
		
		//Check if the input files , paths to ref and sam ok!
		checkParameters(samHybrid, refFile, bamFile);
		
		//create directory for the glf files
		String glfDir = this.OUTPUT_PATH + "glf";
		createDirectory(glfDir);
		
		outFile = glfDir+"/" + IO.replaceFileExtension(bamFile.substring(bamFile.lastIndexOf("/")+1, 
						bamFile.length()), ".glf");
		
	    ArrayList<String> command = new ArrayList<String>();
	    command.add(samHybrid);
	    command.add("pileup");
	    command.add("-g "+ bamFile);
	    command.add("-f " + refFile);
	    command.add("> " + 	outFile);
	    	
	    System.out.println(StringUtils.join(command, " "));
	    Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,logger);
	     		
		
		/**
    	 * OUTPUT
    	 */
     	//Table1
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1_TABLE1, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(outFile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable1 = cont.getTable();


        return new BufferedDataTable[]{outTable1};
		        
    }

    private void checkParameters(String samHybrid, String refFile, String bamFile) throws Exception {
		// check if path is null
        if(bamFile == ""){
        	throw new Exception("No bam file available, something went wrong with the previous node!");
        }
        
        // check path to bam file
        if(!Files.exists(Paths.get(bamFile))){
        	throw new Exception("Path to input bam file: "+bamFile+" does not exist");
        }
        
        String fileextension = PathProcessor.getExt(bamFile);
       
        // check bam format
        if(!fileextension.equals("bam")){
        	throw new Exception("Input file is not in bam format!");
        }
        
        // path to sam hybrid should not be null
        if(samHybrid.equals("")){
        	throw new Exception("Please provide the path to samtools-hybrid");
        }
        
        // check path to samtools-hybrid path
        if(!Files.exists(Paths.get(samHybrid))){
        	throw new Exception("Samtools-hybrid binary: "+samHybrid+" does not exist");
        }
                
        // path to reffile should not be null
        if(refFile.equals("")){
        	throw new Exception("Please provide a reference file");
        }
        
        // check path to reference path
        if(!Files.exists(Paths.get(refFile))){
        	throw new Exception("Reference sequence file: "+refFile+" does not exist");
        }
	}
    
    private void createDirectory(String directory) {
		File dir = new File(directory);

		if (dir.exists())
			dir.delete();
		// if the directory does not exist, create it
		else {
			System.out.println("creating directory: " + dir.toString());
			if (dir.mkdir()) {    
				System.out.println(directory.substring(directory.lastIndexOf("/")+1, directory.length()) + " directory created");  
			}
		}
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

        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {

        // TODO save user settings to the config object.
        
        m_SAMTOOLS_HYBRID.saveSettingsTo(settings);
        m_REF_GENOME.saveSettingsTo(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
        // TODO load (valid) settings from the config object.
        // It can be safely assumed that the settings are valided by the 
        // method below.
        
        m_SAMTOOLS_HYBRID.loadSettingsFrom(settings);
        m_REF_GENOME.loadSettingsFrom(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
        // TODO check if the settings could be applied to our model
        // e.g. if the count is in a certain range (which is ensured by the
        // SettingsModel).
        // Do not actually set any values of any member variables.

        m_SAMTOOLS_HYBRID.validateSettings(settings);
        m_REF_GENOME.validateSettings(settings);

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        
        // TODO load internal data. 
        // Everything handed to output ports is loaded automatically (data
        // returned by the execute method, models loaded in loadModelContent,
        // and user settings set through loadSettingsFrom - is all taken care 
        // of). Load here only the other internals that need to be restored
        // (e.g. data used by the views).

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
       
        // TODO save internal models. 
        // Everything written to output ports is saved automatically (data
        // returned by the execute method, models saved in the saveModelContent,
        // and user settings saved through saveSettingsTo - is all taken care 
        // of). Save here only the other internals that need to be preserved
        // (e.g. data used by the views).

    }

}


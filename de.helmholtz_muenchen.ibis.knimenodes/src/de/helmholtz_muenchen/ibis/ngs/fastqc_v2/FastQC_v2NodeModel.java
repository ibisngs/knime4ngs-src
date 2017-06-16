package de.helmholtz_muenchen.ibis.ngs.fastqc_v2;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.ngs.gatkbaserecalibration.GATKBaseRecalibrationNodeModel;
import de.helmholtz_muenchen.ibis.ngs.vep.VEPNodeModel;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;


/**
 * This is the model implementation of FastQC_v2.
 * 
 *
 * @author Paul Hager
 */
public class FastQC_v2NodeModel extends HTExecutorNodeModel {
	
	// the logger instance
    protected static final NodeLogger LOGGER = NodeLogger.getLogger(FastQC_v2NodeModel.class);
    
    // Node specific params
    public static final String CFGKEY_FASTQC = "Path2FastQC";
    public static final String CFGKEY_OUTFOLDER = "Path2OutFolder";
	public static final String CFGKEY_THREADS = "NumberThreads";
	
	
	private final SettingsModelString m_fastqc = 
			new SettingsModelString(CFGKEY_FASTQC,"");
	
	private final SettingsModelString m_outfolder = 
			new SettingsModelString(CFGKEY_OUTFOLDER,"");
	
	private final SettingsModelIntegerBounded m_threads = 
			new SettingsModelIntegerBounded(CFGKEY_THREADS, 4, 1, Integer.MAX_VALUE);
	
	//The Output Col Names
	public static final String OUT_READFILE1 = "Path2ReadFile1";
	public static final String OUT_READFILE2 = "Path2ReadFile2";
	public static final String OUT_FILTERFILE = "Path2filterfile";
	
	//ReadType: paired-end or single-end
	private static String readType = "";
	
	//IBISKNIMENodesPlugin.FASTQC

	/**
     * Constructor for the node model.
     */
    protected FastQC_v2NodeModel() {
    
        super(1, 1);
        
        addSetting(m_fastqc);
    	addSetting(m_outfolder);
    	addSetting(m_threads);
    	
    	addPrefPageSetting(m_fastqc, IBISKNIMENodesPlugin.FASTQC);
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
    	
    	// Create outfile
    	String outFile1 = this.createOutputName(inFile1);
    	
    	File lockFile = new File(inFile1.substring(0,inFile1.lastIndexOf(".")) + ".FastQC" + SuccessfulRunChecker.LOCK_ENDING);
    	
    	ArrayList<String> cmd = new ArrayList<String>();
    	cmd.add(IO.processFilePath(m_fastqc.getStringValue()));
    	
    	if(!m_outfolder.getStringValue().equals("")){
    		cmd.add("-o="+IO.processFilePath(m_outfolder.getStringValue()));
    	}
    	
    	if(m_threads.getIntValue() != 1){
    		cmd.add("-t="+m_threads.getIntValue());
    	}
    	
    	cmd.add(inFile1);
    	
    	// Parse options into command string array - ensures spaces in path will be correctly interpreted
    	String[] command = new String[cmd.size()];
    	for(int indx=0; indx < cmd.size(); indx++){
    		command[indx] = cmd.get(indx);
    	}
    	
    	LOGGER.info("CMD: "+command.toString());
    	
    	LOGGER.info("Running FastQC...");
		LOGGER.info("Log files can be found in "+inFile1+".stdOut and "+inFile1+".stdErr");
		super.executeCommand(command, outFile1, exec, lockFile, inFile1+".stdOut", inFile1+".stdErr");

		if(readType.equals("paired-end")){
			inFile2 = inData[0].iterator().next().getCell(1).toString();
			
			lockFile = new File(inFile2.substring(0,inFile2.lastIndexOf(".")) + ".FastQC" + SuccessfulRunChecker.LOCK_ENDING);
			command[command.length-1] = inFile2;
			
			LOGGER.info("CMD: "+command.toString());
	    	
	    	LOGGER.info("Running FastQC...");
			LOGGER.info("Log files can be found in "+inFile2+".stdOut and "+inFile2+".stdErr");
			super.executeCommand(command, outFile1, exec, lockFile, inFile2+".stdOut", inFile2+".stdErr");
		}
		
		BufferedDataContainer cont;
	    FileCell[] c;
		
		cont = exec.createDataContainer(createSpecs());
	    
	    
	    if(readType.equals("single-end")){
	    	c = new FileCell[]{
		    		FileCellFactory.create(inFile1)};
	    } else {
	    	c = new FileCell[]{
	    			FileCellFactory.create(inFile1),
	       			FileCellFactory.create(inFile2)};
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
    	if(IO.hasGZipExtension(file)){
    		file = file.replace("\\.fq.gz$", "_trimmed.fq.gz");
    	} else {
    		file = file.replace("\\.fq$", "_trimmed.fq");
    	}
    	
    	String fileName = file.substring(file.lastIndexOf(File.separator), file.length());
    	file = m_outfolder.getStringValue()+fileName;
    	
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
            				new DataColumnSpecCreator(OUT_READFILE1, FastQCell.TYPE).createSpec()});
    	} else {
    		out = new DataTableSpec(
            		new DataColumnSpec[]{
            				new DataColumnSpecCreator(OUT_READFILE1, FastQCell.TYPE).createSpec(),
            				new DataColumnSpecCreator(OUT_READFILE2, FastQCell.TYPE).createSpec()});
    	}
    	
    	return out;
    }
    
    /**
     * @see org.knime.core.node.NodeModel
     *      #validateSettings(org.knime.core.node.NodeSettingsRO)
     */
     @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
    	// delegate this to the settings models
    	
    	 m_fastqc.validateSettings(settings);
    	 m_outfolder.validateSettings(settings);
    	 m_threads.validateSettings(settings);
    }
     
     /**
      * @see org.knime.core.node.NodeModel
      *      #loadValidatedSettingsFrom(org.knime.core.node.NodeSettingsRO)
      */
      @Override
     protected void loadValidatedSettingsFrom(
             final NodeSettingsRO settings)
             throws InvalidSettingsException {
             
     	// loads the values from the settings into the models.
         // It can be safely assumed that the settings are validated by the 
         // method below.
         
    	  m_fastqc.loadSettingsFrom(settings);
    	  m_outfolder.loadSettingsFrom(settings);
    	  m_threads.loadSettingsFrom(settings);

     }
      
      /**
       * @see org.knime.core.node.NodeModel
       *      #saveSettingsTo(org.knime.core.node.NodeSettings)
       */
       @Override
      protected void saveSettingsTo(final NodeSettingsWO settings) {

          // save settings to the config object.
      	
    	   m_fastqc.saveSettingsTo(settings);
    	   m_outfolder.saveSettingsTo(settings);
    	   m_threads.saveSettingsTo(settings);
      }
 
}


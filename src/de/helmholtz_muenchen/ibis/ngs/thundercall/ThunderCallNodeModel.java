package de.helmholtz_muenchen.ibis.ngs.thundercall;

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
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.lofs.PathProcessor;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;


/**
 * This is the model implementation of ThunderCall.
 * 
 *
 * @author Tanzeem Haque
 */
public class ThunderCallNodeModel extends NodeModel {
    
    // the logger instance
    private static final NodeLogger logger = NodeLogger.getLogger(ThunderCallNodeModel.class);
        
    /**
	 * Config Keys
	 */
	public static final String CFGKEY_THUNDER_PATH = "Thunder path";
	public static final String CFGKEY_BASE_NAME = "Base name";
	public static final String CFGKEY_POST_PROB ="Posterior probabilty";
	public static final String CFGKEY_MIN_DEPTH = "Mininum depth";
	public static final String CFGKEY_MAX_DEPTH = "Maximum depth";

	private String OUTPUT_PATH;
	
	/**
	 * Initial default values
	 */
	static final double DEFAULT_POST_PROB = 0.9;
	static final int DEFAULT_MIN_DEPTH = 10;
	static final int DEFAULT_MAX_DEPTH = 1000;

	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_THUNDER = new SettingsModelString(CFGKEY_THUNDER_PATH, "");
    private final SettingsModelString m_BASE_NAME = new SettingsModelString(CFGKEY_BASE_NAME, "");
    private final SettingsModelDoubleBounded m_POST_PROB = new SettingsModelDoubleBounded(CFGKEY_POST_PROB, DEFAULT_POST_PROB, 0.0, 1.0);
    private final SettingsModelIntegerBounded m_MIN_DEPTH = new SettingsModelIntegerBounded(CFGKEY_MIN_DEPTH, DEFAULT_MIN_DEPTH, 1, 100);
    private final SettingsModelIntegerBounded m_MAX_DEPTH = new SettingsModelIntegerBounded(CFGKEY_POST_PROB, DEFAULT_MAX_DEPTH, 100, 10000);

  //The Output Col Names
  	public static final String OUT_COL1_TABLE1 = "OUTFILE";
    /**
     * Constructor for the node model.
     */
    protected ThunderCallNodeModel() {
    
        // TODO one incoming port and one outgoing port is assumed
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

        // TODO do something here
        
        String thunder = m_THUNDER.getStringValue();
        String baseName = m_BASE_NAME.getStringValue();
        double postProb = m_POST_PROB.getDoubleValue();
        int minDepth = m_MIN_DEPTH.getIntValue();
        int maxDepth = m_MAX_DEPTH.getIntValue();
        String outFile = "";

        /**
         * Connect with parallel chunk and set the number of chunks accordingly
         * Will always work with a trio
         * chr21_m
         * chr21_f
         * chr21_c
         * 
         * then the next chunk would be
         * chr22_m
         * chr22_f
         * chr22_c
         */
    	ArrayList<String> glf_arrList = new ArrayList<>(3);

        CloseableRowIterator it = inData[0].iterator();
		while (it.hasNext()) {
			DataRow row = it.next();
			glf_arrList.add(row.getCell(0).toString());
		}
		
		this.OUTPUT_PATH = glf_arrList.get(0).substring(0, glf_arrList.get(0).lastIndexOf("/"))+"/";

		// path to thunder should not be null
        if(thunder.equals("")){
        	throw new Exception("Please provide the path to thunder_GPT_Freq tool");
        }
        
        // check path to thunder path
        if(!Files.exists(Paths.get(thunder))){
        	throw new Exception("Thunder_GPT_Freq binary: "+thunder+" does not exist");
        }        
		
		//create directoty for the Thunder Outputs
		String thunderDir = this.OUTPUT_PATH + "ThunderCall";
		String glfDir = thunderDir.replace("ThunderCall", "glf");
		createDirectory(thunderDir);
		/**
		 * /home/ibis/tanzeem.haque/Documents/3rd_party_tools/thunder/GPT_Freq_V011/GPT_Freq 
		 * -b /home/ibis/tanzeem.haque/Documents/3rd_party_tools/thunder/thunder_011/examples/testTZ/tin/Q12 
		 * -p 0.9 
		 * --minDepth 1 
		 * --maxDepth 30 glf/*.glf 
		 * --stats
		 * > /home/ibis/tanzeem.haque/Documents/3rd_party_tools/thunder/thunder_011/examples/testTZ/tin/GPT.Q12.log
		 */

		String chr = glf_arrList.get(0).split("_")[0];
		ArrayList<String> command = new ArrayList<String>();
    	command.add(thunder);
    	command.add("-b "+ thunderDir + "/" + baseName);
    	command.add("-p " + postProb);
    	command.add("--minDepth " + minDepth);
    	command.add("--maxDepth " + maxDepth);
    	command.add(glfDir + "/" + chr + "*.glf");
    	command.add("> " + 	thunderDir + "/" + "GPT."+baseName+".log");
    	
    	System.out.println(StringUtils.join(command, " "));
    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,logger);
     	
    	outFile = thunderDir + "/" + baseName + "." + chr;
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

	/**
	 * 
	 */
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
        m_THUNDER.saveSettingsTo(settings);
        m_POST_PROB.saveSettingsTo(settings);
        m_BASE_NAME.saveSettingsTo(settings);
        m_MIN_DEPTH.saveSettingsTo(settings);
        m_MAX_DEPTH.saveSettingsTo(settings);

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
        
        m_THUNDER.loadSettingsFrom(settings);
        m_POST_PROB.loadSettingsFrom(settings);
        m_BASE_NAME.loadSettingsFrom(settings);
        m_MIN_DEPTH.loadSettingsFrom(settings);
        m_MAX_DEPTH.loadSettingsFrom(settings);


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

        m_THUNDER.validateSettings(settings);
        m_POST_PROB.validateSettings(settings);
        m_BASE_NAME.validateSettings(settings);
        m_MIN_DEPTH.validateSettings(settings);
        m_MAX_DEPTH.validateSettings(settings);


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


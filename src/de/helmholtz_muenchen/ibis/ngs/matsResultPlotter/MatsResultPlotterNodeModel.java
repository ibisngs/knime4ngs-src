package de.helmholtz_muenchen.ibis.ngs.matsResultPlotter;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.ngs.matsResultIndexer.MatsResultIndexerNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode.ExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of MatsResultPlotter.
 * 
 *
 * @author Michael Kluge
 */
public class MatsResultPlotterNodeModel extends ExecutorNodeModel {
	
	public static final String OUT_COL1 = "OutputFolder";
	public static final String PYTHON_BIN = "python";
	
	 // keys for SettingsModels
	protected static final String CFGKEY_BINARY_PATH 		= "BinaryPathSashimi";
    protected static final String CFGKEY_OUTPUT_FILE 		= "OutputFolder";
    protected static final String CFGKEY_INPUT_FILE			= "InputFolder";
    protected static final String CFGKEY_SETTINGS_FILE		= "SettingsFolder";
    protected static final String CFGKEY_FDR				= "FDRFilter";
    protected static final String CFGKEY_INCLUDE_READS		= "IncludeReads";
    
    // initial default values for SettingsModels
    protected static final String DEFAULT_BINARY_PATH 		= "";
    protected static final String DEFAULT_OUTPUT_FOLDER 	= "./output/";	
    protected static final String DEFAULT_INPUT_FOLDER 		= "";
    protected static final String DEFAULT_SETTINGS_FILE		= "";
    protected static final double DEFAULT_FDR				= 0.01;
    protected static final boolean DEFAULT_INCLUDE_READS	= true;
    
    private final static String NAME_OF_PLOT_EVENT 			= "--plot-event";
    private final static String NAME_OF_OUTPUT_FOLDER 		= "--output-dir";	
    
    private final SettingsModelString SET_BINARY_PATH	= new SettingsModelString(CFGKEY_BINARY_PATH, DEFAULT_BINARY_PATH);
    private final SettingsModelString SET_OUTPUT_FILE = new SettingsModelString(CFGKEY_OUTPUT_FILE, DEFAULT_OUTPUT_FOLDER);
    private final SettingsModelString SET_INPUT_FILE = new SettingsModelString(CFGKEY_INPUT_FILE, DEFAULT_INPUT_FOLDER);
    private final SettingsModelString SET_SETTINGS_FILE = new SettingsModelString(CFGKEY_SETTINGS_FILE, DEFAULT_SETTINGS_FILE);
    private final SettingsModelDouble SET_FDR = new SettingsModelDouble(CFGKEY_FDR, DEFAULT_FDR);
    private final SettingsModelBoolean SET_INCLUDE_READS = new SettingsModelBoolean(CFGKEY_INCLUDE_READS, DEFAULT_INCLUDE_READS);
    
    // the logger instance
	private static final NodeLogger LOGGER = NodeLogger.getLogger(MatsResultIndexerNodeModel.class);
	
	/**
	 * Constructor
	 */
	protected MatsResultPlotterNodeModel() {
		super(1, 1, true, true);
		init();
	}
	 
	 /**
     * {@inheritDoc}
     */
    @Override
    public void init() {    	
    	this.addSetting(SET_INPUT_FILE);
    	this.addSetting(SET_BINARY_PATH);
    	this.addSetting(SET_OUTPUT_FILE);
    	this.addSetting(SET_SETTINGS_FILE);
    	this.addSetting(SET_FDR);
    	this.addSetting(SET_INCLUDE_READS);
    }
    
	@Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
		return new DataTableSpec[]{this.getDataOutSpec1()};
	}
	
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
		// get input data
		String indexFolder = inData[0].iterator().next().getCell(0).toString();
		String inputFolder = SET_INPUT_FILE.getStringValue();
		String outputFile = SET_OUTPUT_FILE.getStringValue();
		String binaryfile = SET_BINARY_PATH.getStringValue();
		String settingsFile = SET_SETTINGS_FILE.getStringValue();
		double FDR = SET_FDR.getDoubleValue();
		boolean includeReads = SET_INCLUDE_READS.getBooleanValue();
		
		// check, if settings are valid
		if(FDR > 1 || FDR < 0)
			throw new IllegalArgumentException("FDR must be between 0 and 1.");
		if(!new File(indexFolder).exists())
			throw new IllegalArgumentException("Index folder '" + indexFolder + "' was not found");
		if(!new File(settingsFile).exists())
			throw new IllegalArgumentException("Settings file '" + settingsFile + "' was not found");
		if(!new File(binaryfile).exists())
			throw new IllegalArgumentException("Sashimi plot script '" + settingsFile + "' was not found");
		
		// create folder, if not existent
		File out = new File(outputFile);
		if(!out.exists())
		 out.mkdirs();
		
		// find events that must be ploted
		if(new File(inputFolder).exists() & new File(inputFolder).isDirectory()) {
			FDRFinder finder = new FDRFinder(inputFolder, includeReads, FDR);
			HashMap<String, ArrayList<String>> events = finder.getSigEvents();
			
			for(String type : events.keySet()) {
				ArrayList<String> eve = events.get(type);
				System.out.println("type; " + type);
				for(String e : eve) {
					System.out.println(e);
				}
			}
			
			// build command array
			ArrayList<String> command = new ArrayList<String>();
			command.add(PYTHON_BIN);
			command.add(binaryfile);
			command.add(NAME_OF_PLOT_EVENT);
			command.add(null);
			command.add(indexFolder);
			command.add(settingsFile);
			command.add(NAME_OF_OUTPUT_FOLDER);
			command.add(null);
			
			String[] c = command.toArray(new String[0]);
			
			for(String type : events.keySet()) {
				// create fold for event
				File outSub = new File(out.getAbsolutePath() + File.separator + type);
				if(!outSub.exists()) 
					outSub.mkdir();
				c[7] = outSub.getAbsolutePath();
				// plot the events
				for(String name : events.get(type)) {
					LOGGER.info("plotting: " + type + "  " + name + "...");
					c[3] = name;
					
					// execute the command.
					String stdOutFile = c[7] + File.separator + name + ".log";
					Executor.executeCommand(c, exec, LOGGER, stdOutFile);
				}
			}
			 return getOutputData(exec, "multiple commands, see log file", inData); 
		}
		else {
			throw new IllegalArgumentException("Folder '"+ inputFolder +"' was not found.");
		}
	}
	 
    /**
     * returns the first output specifications of this node
     * @return
     */
    private DataTableSpec getDataOutSpec1() {
    	return new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, StringCell.TYPE).createSpec()});
    }

    /**
     * gets the output data
     * @param exec
     * @param command
     * @param inData
     * @return
     */
	protected BufferedDataTable[] getOutputData(ExecutionContext exec, String command, BufferedDataTable[] inData) {
		BufferedDataContainer cont = exec.createDataContainer(getDataOutSpec1());
		
    	DataCell[] c = new DataCell[]{
    			new StringCell(SET_OUTPUT_FILE.getStringValue())};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
        return new BufferedDataTable[]{cont.getTable()};
	}
}


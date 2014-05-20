package de.helmholtz_muenchen.ibis.utils.abstractNodes.StatisticMerger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SettingsStorageNodeModel;

/**
 * Default implementation of a "Statistic Merger Node"
 * @author Michael Kluge
 *
 */
public abstract class StatisticMergerNodeModel extends SettingsStorageNodeModel {
	
	// Available module names
	private static final String MODULE_ALL = "Merge ALL Modules";
	public static final String FILE_ENDING = ".txt";
	protected static final LinkedHashMap<String, ArrayList<String>> MODULE_NAMES = new LinkedHashMap<String, ArrayList<String>>();
	
	 // keys for SettingsModels
	protected static final String CFGKEY_INPUT_FOLDER 	= "InputFolder";
	protected static final String CFGKEY_OUTPUT_FOLDER 	= "OutputFolder";
	protected static final String CFGKEY_MODULE_ALL 	=  MODULE_ALL;
	
    // initial default values for SettingsModels    
	protected static final String DEFAULT_INPUT_FOLDER 	= "";		
	protected static final String DEFAULT_OUTPUT_FOLDER 	= "";	
    
    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_INPUT_FOLDER	= new SettingsModelString(CFGKEY_INPUT_FOLDER, DEFAULT_INPUT_FOLDER);
    private final SettingsModelString SET_OUTPUT_FOLDER = new SettingsModelString(CFGKEY_OUTPUT_FOLDER, DEFAULT_OUTPUT_FOLDER);
    private final SettingsModelBoolean SET_ALL_MODULES = new SettingsModelBoolean(CFGKEY_MODULE_ALL, false);
	
    /**
     * Constructor for the node model.
     */
    protected StatisticMergerNodeModel() {
        super(0, 2);
        init();
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public void init() {
    	addSetting(SET_INPUT_FOLDER);
    	addSetting(SET_OUTPUT_FOLDER);
    	
    	// add only, if more than one module is there
    	if(this.getModuleNames().size() > 1)
    		addSetting(SET_ALL_MODULES);
    	
        // add modules
        for(String moduleName : this.getModuleNames()) {
        	addSetting(new SettingsModelBoolean(moduleName, false));
        }
    }
    

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
    	// create output folder
    	File outputFolder = new File(SET_OUTPUT_FOLDER.getStringValue());
    	if(!outputFolder.isDirectory())
    			outputFolder.mkdirs();
    	
    	// test if input folder is there
    	String inputPath = SET_INPUT_FOLDER.getStringValue();
    	if(!new File(inputPath).isDirectory()) {
    		throw new IllegalArgumentException("Input path '" + inputPath + "' was not found.");
    	}
    		
    	// get all files in folder
    	ArrayList<String> dataFiles = this.getStatisticFiles(inputPath);
    	BufferedDataContainer cont1 = exec.createDataContainer(getDataOutSpec1());
    	BufferedDataContainer cont2 = exec.createDataContainer(getDataOutSpec2());
    	int i = 0;
    	
    	// run through all activated modules
    	for(String moduleName : this.getModuleNames()) {
    		SettingsModelBoolean option = (SettingsModelBoolean) getSetting(moduleName);
    		if(!moduleName.equals(MODULE_ALL)) {
	    		// check, if that module was activated in the GUI
	    		if(option.getBooleanValue()) {
			    	String moduleFileName = moduleName.replace(" ", "_") + FILE_ENDING;
			    	boolean writeHeader = true;
			    	String outfile = new File(SET_OUTPUT_FOLDER.getStringValue() + File.separator + moduleFileName).getAbsolutePath();
			    	BufferedWriter outfileBW = new BufferedWriter(new FileWriter(outfile));
			    	
			    	// run though all files
			    	for(String file : dataFiles) {
			    		extractTable(moduleName, file, outfileBW, writeHeader);
			    		writeHeader = false;
			    	}
			    	
			    	outfileBW.flush();
			    	finalize(outfileBW);
			    	// close outfile
			    	outfileBW.close();
			    	// write output Table
			    	DataCell[] c = new DataCell[]{
			    			new StringCell(moduleName),
			    			new StringCell(outfile)};
			    	cont1.addRowToTable(new DefaultRow("Row"+i, c));
			    	i++;
	    		}
    		}
    	}
    	cont1.close();
    	
    	// write second table
    	i = 0;
    	for(String file : dataFiles) {
    		DataCell[] c = new DataCell[]{
	    			new StringCell(file)};
	    	cont2.addRowToTable(new DefaultRow("Row"+i, c));
	    	i++;
    	}
    	cont2.close();
		
        return new BufferedDataTable[]{cont1.getTable(), cont2.getTable()};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {}

    /**
     * returns the first output specifications of this node
     * @return
     */
    private DataTableSpec getDataOutSpec1() {
    	return new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator("Module", StringCell.TYPE).createSpec(),
    					new DataColumnSpecCreator("OutputPath", StringCell.TYPE).createSpec()});
    }
    
    /**
     * returns the second output specifications of this node
     * @return
     */
    private DataTableSpec getDataOutSpec2() {
    	return new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator("InputFile", StringCell.TYPE).createSpec()});
    }
    
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	// test if input folder is set
    	if(SET_INPUT_FOLDER.getStringValue().isEmpty())
    		throw new InvalidSettingsException("Input folder path may not be empty.");    	
    	// test if output folder is set
    	if(SET_OUTPUT_FOLDER.getStringValue().isEmpty())
		throw new InvalidSettingsException("Output folder path may not be empty.");
    	
    	boolean activated = false;
    	// test if at least one module is activated
    	for(String moduleName : this.getModuleNames()) {
    		SettingsModelBoolean option = (SettingsModelBoolean) getSetting(moduleName);
	    	if(option.getBooleanValue()) {
	    		activated = true;
	    		break;
	    	}
    	}
    	if(!activated)
    		throw new InvalidSettingsException("At least one module must be activated.");
    	
        return new DataTableSpec[]{getDataOutSpec1(), getDataOutSpec2()};
    }

	@Override
	protected void loadInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {
	}

	@Override
	protected void saveInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {
	}
	
	/**
	 * Gets all statistic files in a folder and in subfolders
	 * @param path Path to start with
	 * @return
	 */
	protected ArrayList<String> getStatisticFiles(String path) {
		return getStatisticFiles(path, new ArrayList<String>());
	}
	
	
	/**
	 * Gets all statistic files in a folder and in subfolders
	 * @param path Path to start with
	 * @param files Files which were already found
	 * @return
	 */
	private ArrayList<String> getStatisticFiles(String path, ArrayList<String> files) {
		File fpath = new File(path);

		// test all files in folder
		for(File f : fpath.listFiles()) {
			// if name is ok, add it
			if(f.isFile() && f.getName().matches(this.getNameOfStatisticFile())) {
				files.add(f.getAbsolutePath());
			}
			// call function recursive
			else if(f.isDirectory()) {
				getStatisticFiles(f.getAbsolutePath(), files);
			}
		}
		return files;
	}
	
	/**
	 * adds a new module names
	 * @param mergerName
	 * @param moduleNames
	 */
	public static void addModuleNames(String mergerName, ArrayList<String> moduleNames) {
		MODULE_NAMES.put(mergerName, moduleNames);
	}

	/**
	 * get the module names for a specific merger
	 * @param mergerName
	 * @return
	 */
    public static ArrayList<String> getModuleNames(String mergerName) {
    	return MODULE_NAMES.get(mergerName);
    }
    
    /**
     * Names of the modules
     * @return
     */
    public ArrayList<String> getModuleNames() {
    	return getModuleNames(this.getMergerName());
    }
	
    /****************************** ABSTRACT METHODS **********************************/
    
	/**
	 * Extracts the needed data out of the file for a specific module
	 * @param moduleName name of the module
	 * @param file Input file
	 * @param outfile output file
	 * @param writeHeader true, if header must be written
	 * @return
	 */
	protected abstract boolean extractTable(String moduleName, String file, BufferedWriter outfile, boolean writeHeader);
	
    
    /**
     * Name of the file where the statistic should be extracted
     * @return
     */
    public abstract String getNameOfStatisticFile();
    
    /**
     * Name of the statistic merger
     * @return
     */
    public abstract String getMergerName();
    
    /**
     * Is called after all files for a module were processed
     * @param outfile
     */
	public abstract void finalize(BufferedWriter outfile) throws IOException;
}
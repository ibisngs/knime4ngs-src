package de.helmholtz_muenchen.ibis.ngs.fastqcStatisticMerger;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

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
 * This is the model implementation of FastqcStatisticMergerNodeModel.
 * This node can be used to merge the statistics of the FASTQC Node
 *
 * @author Michael Kluge
 */
public class FastqcStatisticMergerNodeModel extends SettingsStorageNodeModel {
    
	public static final String DATA_FILE = "fastqc_data.txt";
	public static final String MODULE_PREFIX = ">>";
	public static final String MODULE_END = "END_MODULE";
	public static final String BASIC_FILENAME = "Filename";
	public static final String BASIC_MODULE = "Basic Statistics";
	public static final String TAB = "\t";
	public static final String HEADER_PREFIX = "#";
	public static final String HEADER_NAME = "Filename";
	public static final String HEADER_STATUS = "#Module" + TAB + "Status";
	public static final String FILE_ENDING = ".txt";
	protected static final String KEY = "MODULE_KEY_";
	
	// Available module names
	private static final String MODULE_BASIC = "Basic Statistics";
	private static final String MODULE_BSQ = "Per base sequence quality";
	private static final String MODULE_SQS = "Per sequence quality scores";
	private static final String MODULE_BSC = "Per base sequence content";
	private static final String MODULE_BGCC = "Per base GC content";
	private static final String MODULE_SGCC = "Per sequence GC content";
	private static final String MODULE_BNC = "Per base N content";
	private static final String MODULE_SLD = "Sequence Length Distribution";
	private static final String MODULE_SDL = "Sequence Duplication Levels";
	private static final String MODULE_OS = "Overrepresented sequences";
	private static final String MODULE_FOC = "Filter-Options Collector";
	private static final String MODULE_STATUS = "Merge Status of all Modules";
	private static final String MODULE_ALL = "Merge ALL Modules";
	protected static final ArrayList<String> MODULE_NAMES = new ArrayList<String>();
	
	 // keys for SettingsModels
	protected static final String CFGKEY_INPUT_FOLDER 	= "InputFolder";
	protected static final String CFGKEY_OUTPUT_FOLDER 	= "OutputFolder";
	protected static final String CFGKEY_MODULE_ALL 	=  MODULE_ALL;
	
    // initial default values for SettingsModels    
	private static final String DEFAULT_INPUT_FOLDER 	= "";		
    private static final String DEFAULT_OUTPUT_FOLDER 	= "";	
    
    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_INPUT_FOLDER	= getSettingsModelString(CFGKEY_INPUT_FOLDER, this);
    private final SettingsModelString SET_OUTPUT_FOLDER = getSettingsModelString(CFGKEY_OUTPUT_FOLDER, this);

	 /**
     * add the used settings
     */
    static {
        // add values for SettingsModelString
        addSettingsModelString(CFGKEY_INPUT_FOLDER, DEFAULT_INPUT_FOLDER);
        addSettingsModelString(CFGKEY_OUTPUT_FOLDER, DEFAULT_OUTPUT_FOLDER);
        
        // add modules
        MODULE_NAMES.add(MODULE_BASIC);
        MODULE_NAMES.add(MODULE_BSQ);
        MODULE_NAMES.add(MODULE_SQS);
        MODULE_NAMES.add(MODULE_BSC);
        MODULE_NAMES.add(MODULE_BGCC);
        MODULE_NAMES.add(MODULE_SGCC);
        MODULE_NAMES.add(MODULE_BNC);
        MODULE_NAMES.add(MODULE_SLD);
        MODULE_NAMES.add(MODULE_SDL);
        MODULE_NAMES.add(MODULE_OS);
        MODULE_NAMES.add(MODULE_FOC);
        MODULE_NAMES.add(MODULE_STATUS);
        MODULE_NAMES.add(MODULE_ALL);
        
        // add the SettingsModelString
        for(String moduleName : MODULE_NAMES) {
        	addSettingsModelBoolean(moduleName, false);
        }
    }
	
	
    /**
     * Constructor for the node model.
     */
    protected FastqcStatisticMergerNodeModel() {
        super(0, 2);
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

    	// get all files in folder
    	ArrayList<String> dataFiles = getFastqcSummaryFiles(SET_INPUT_FOLDER.getStringValue());
    	BufferedDataContainer cont1 = exec.createDataContainer(getDataOutSpec1());
    	BufferedDataContainer cont2 = exec.createDataContainer(getDataOutSpec2());
    	int i = 0;
    	
    	// run through all activated modules
    	for(String moduleName : MODULE_NAMES) {
    		SettingsModelBoolean option = getSettingsModelBoolean(moduleName);
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
    	for(String moduleName : MODULE_NAMES) {
    		SettingsModelBoolean option = getSettingsModelBoolean(moduleName);
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
	 * Gets all fastq data files in a folder and in subfolders
	 * @param path Path to start with
	 * @return
	 */
	protected ArrayList<String> getFastqcSummaryFiles(String path) {
		return getFastqcSummaryFiles(path, new ArrayList<String>());
	}
	
	
	/**
	 * Gets all fastq data files in a folder and in subfolders
	 * @param path Path to start with
	 * @param files Files which were already found
	 * @return
	 */
	private ArrayList<String> getFastqcSummaryFiles(String path, ArrayList<String> files) {
		File fpath = new File(path);

		// test all files in folder
		for(File f : fpath.listFiles()) {
			// if name is ok, add it
			if(f.isFile() && f.getName().equals(DATA_FILE)) {
				files.add(f.getAbsolutePath());
			}
			// call function recursive
			else if(f.isDirectory()) {
				getFastqcSummaryFiles(f.getAbsolutePath(), files);
			}
		}
		return files;
	}
	
	/**
	 * Extracts the needed data out of the file for a specific module
	 * @param moduleName name of the module
	 * @param file Input file
	 * @param outfile output file
	 * @param writeHeader true, if header must be written
	 * @return
	 */
	protected boolean extractTable(String moduleName, String file, BufferedWriter outfile, boolean writeHeader) {
		boolean isStatusMode = MODULE_STATUS.equals(moduleName);
		File f = new File(file);
		// test if correct file is given
		if(!f.isFile() || !f.canRead() || !f.getName().equals(DATA_FILE))
			return false;
		
		// open file
		try {
			String line;
			BufferedReader r = new BufferedReader(new FileReader(f));

			String name = null;
			boolean hasFoundBasicModule = false;
			boolean hasFoundCorrectModule = false;

			// read lines
			while((line = r.readLine()) != null) {
				// find start of basic module
				if(!hasFoundBasicModule && !hasFoundCorrectModule) {
					if(line.startsWith(MODULE_PREFIX + BASIC_MODULE)) 
						hasFoundBasicModule = true;
				}
				// get filename of fasta file which was analysed with FASTQ
				else if(hasFoundBasicModule && !hasFoundCorrectModule) {
					// find name of file
					if(name == null && line.startsWith(BASIC_FILENAME)) {
						name = line.replace(BASIC_FILENAME, "").replace(TAB, "");
						// go back to start of file (if basic module should be merged)
						r = new BufferedReader(new FileReader(f));
						// check, if status should be merged
						if(isStatusMode) 
							hasFoundCorrectModule = true; // all modules are correct ;)
					}
					else if(line.startsWith(MODULE_PREFIX + moduleName)) 
						hasFoundCorrectModule = true;
				}
				// is in correct module --> write lines
				else if(hasFoundCorrectModule) {
					if(!isStatusMode && line.equals(MODULE_PREFIX + MODULE_END)) 
						break;
					else {
						// check, if starts with char for header line
						if(line.startsWith(HEADER_PREFIX)) {
							if(writeHeader) {
								if(!isStatusMode)
									outfile.write(line + TAB + HEADER_NAME);
								else
									outfile.write(HEADER_STATUS + TAB + HEADER_NAME);
								// do not write another header in case of status merger
								outfile.newLine();
								writeHeader = false;
							}
						}
						else if(!isStatusMode) {
							outfile.write(line + TAB + name);
							outfile.newLine();
						}
						// statusMode
						else if(line.startsWith(MODULE_PREFIX) && !line.equals(MODULE_PREFIX + MODULE_END)) {
							String currentModuleName = line.replaceFirst(MODULE_PREFIX, "").split(TAB)[0];
							outfile.write(currentModuleName);
							outfile.write(TAB);
							outfile.write(line.replaceFirst(MODULE_PREFIX + currentModuleName + TAB, "")); 
							outfile.write(TAB);
							outfile.write(name);
							outfile.newLine();
						}
					}
				}
			}
			outfile.flush();
			// close the file
			r.close();
			return true;
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		return false;
	}
}


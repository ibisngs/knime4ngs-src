package de.helmholtz_muenchen.ibis.utils.abstractNodes.FileSelector;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;

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
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SettingsStorageNodeModel;


/**
 * <code>FileSelectorNodeModel</code> can be extended to select multiple files of a specific type.
 *
 * @author Michael Kluge
 */
public abstract class FileSelectorNodeModel extends SettingsStorageNodeModel {
	   
    // keys for SettingsModels
    protected static final String CFGKEY_FILE_LIST 			= "FileList";
    protected static final String CFGKEY_FILE_DIR 			= "FileDir";
    protected static final String CFGKEY_FILE_FILE 			= "FileFile";
    protected static final String CFGKEY_FILE_LIST_DISPLAY 	= "FileListDisplay";
    protected static final String CFGKEY_REGEX			 	= "FilenameRegex";
    
    // initial default values for SettingsModels
    protected static final String DEFAULT_FILE_DIR = "-";
    protected static final String DEFAULT_FILE_FILE = "-";
    protected static final String DEFAULT_REGEX = ".*";
    
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_FILE_DIR 					= new SettingsModelString(CFGKEY_FILE_DIR, DEFAULT_FILE_DIR);
    private final SettingsModelString SET_FILE_FILE 				= new SettingsModelString(CFGKEY_FILE_FILE, DEFAULT_FILE_FILE);
    private final SettingsModelString SET_NAME_REGEX				= new SettingsModelString(CFGKEY_REGEX, DEFAULT_REGEX);
    
    // storage 
    private final HashSet<String> FILES		= new HashSet<String>();
    private boolean hasConfigureOpendOnce 	= false; // true, if configure was opend once
    
    // the logger instance
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FileSelectorNodeModel.class);
   	
    /**
     * Constructor for the node model.
     */
    protected FileSelectorNodeModel() {
        super(0, 1);
        addSetting(SET_FILE_DIR);
		addSetting(SET_FILE_FILE);
		addSetting(SET_NAME_REGEX);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	// check, if some files where selected
    	if(hasConfigureOpendOnce && FILES.size() == 0)
    		throw new InvalidSettingsException("Select at least one " + getFiletypeName() + " file.");
    	
    	return new DataTableSpec[]{getDataOutSpec1()};
    }
    
    /**
     * returns the first output specifications of this node
     * @return
     */
    private DataTableSpec getDataOutSpec1() {
    	DataColumnSpecCreator col1 = new DataColumnSpecCreator(getNameOfOutputCol(), StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec()};
        
        return new DataTableSpec(cols);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
    	exec.setProgress(0.01);

    	/******************* PREPARE OUTPUT **********************/
    	DataTableSpec table = getDataOutSpec1();
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	
    	// write output to table
    	int i = 0;
    	for(String filename : this.FILES)
    	{
	    	StringCell cl1 = new StringCell(filename);
	    	DataCell[] c = new DataCell[]{cl1};
	    	DefaultRow r = new DefaultRow("Row"+i, c);
	    	cont.addRowToTable(r);
	    	i++;
    	}
    	
    	cont.close();
    	BufferedDataTable out = cont.getTable();
    	
    	exec.setProgress(1.0); // we are ready! ;)
        return new BufferedDataTable[]{out};
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
    	this.FILES.clear();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
    	// clean the old data
    	FILES.clear();
    	// check, if data is set
        if (settings.containsKey(FileSelectorNodeModel.CFGKEY_FILE_LIST)) {
        	try {
        		// add the values
				for(String s : settings.getStringArray(FileSelectorNodeModel.CFGKEY_FILE_LIST))
					this.FILES.add(s);
			} catch (InvalidSettingsException e) {
				LOGGER.error(e.getStackTrace());
			}
        }
    }
    
    @Override
	protected void saveSettingsTo(final NodeSettingsWO settings) {
    	super.saveSettingsTo(settings);
 
    	settings.addStringArray(FileSelectorNodeModel.CFGKEY_FILE_LIST, FILES.toArray(new String[FILES.size()]));
	}

	@Override
	protected void validateSettings(NodeSettingsRO settings) throws InvalidSettingsException {
		super.validateSettings(settings);
        // configure must have been opened or we won't be here
        hasConfigureOpendOnce = true;
	}
	
	@Override
	protected void loadInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {	
	}

	@Override
	protected void saveInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {
	}
	
	/************************************** ABSTRACT METHODS **********************************************/
	/******************************************************************************************************/
	
    /**
     * Must return the name of the output col for the specific file selector
     * @return
     */
    public abstract String getNameOfOutputCol();
    
    /**
     * Returns the filetype how it is named in the GUI
     * @return
     */
    protected abstract String getFiletypeName();
}


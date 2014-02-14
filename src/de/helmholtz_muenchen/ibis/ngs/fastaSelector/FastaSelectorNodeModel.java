package de.helmholtz_muenchen.ibis.ngs.fastaSelector;

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

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SettingsStorageNodeModel;


/**
 * <code>FastaSelectorNodeModel</code> can be used to select multiple fasta files.
 *
 * @author Michael Kluge
 */
public class FastaSelectorNodeModel extends SettingsStorageNodeModel {
    
    // name of the output variables
    public static String OUTPUT_NAME_FASTA_FILES = "FastaFile";
    
    // keys for SettingsModels
    protected static final String CFGKEY_FASTA_LIST 		= "FastaList";
    protected static final String CFGKEY_FASTA_DIR 			= "FastaDir";
    protected static final String CFGKEY_FASTA_FILE 		= "FastaFile";
    protected static final String CFGKEY_FASTA_LIST_DISPLAY = "FastaListDisplay";

    // initial default values for SettingsModels
    protected static final String DEFAULT_FASTA_DIR = "-";
    protected static final String DEFAULT_FASTA_FILE = "-";
       
    // storage 
    private final HashSet<String> FASTA_FILES = new HashSet<String>();
    
    // the logger instance
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FastaSelectorNodeModel.class);
   
    /**
     * Constructor for the node model.
     */
    protected FastaSelectorNodeModel() {
        super(0, 1);
        
        // add values for SettingsModelString
        addSettingsModelString(CFGKEY_FASTA_DIR, DEFAULT_FASTA_DIR);
        addSettingsModelString(CFGKEY_FASTA_FILE, DEFAULT_FASTA_FILE);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
    	
    	exec.setProgress(0.1);

    	/******************* PREPARE OUTPUT **********************/
		DataColumnSpecCreator col1 = new DataColumnSpecCreator(OUTPUT_NAME_FASTA_FILES, StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec()};
    	DataTableSpec table = new DataTableSpec(cols);
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	
    	// write output to table
    	int i = 0;
    	for(String fastaFile : this.FASTA_FILES)
    	{
	    	StringCell cl1 = new StringCell(fastaFile);
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
    	this.FASTA_FILES.clear();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {

        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	// clean the old data
    	this.FASTA_FILES.clear();
    	// check, if data is set
        if (settings.containsKey(FastaSelectorNodeModel.CFGKEY_FASTA_LIST)) {
        	try {
        		// add the values
				for(String s : settings.getStringArray(FastaSelectorNodeModel.CFGKEY_FASTA_LIST))
					this.FASTA_FILES.add(s);
			} catch (InvalidSettingsException e) {
				LOGGER.error(e.getStackTrace());
			}
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	
    }

	@Override
	protected void loadInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {

	}

	@Override
	protected void saveInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {

	}
}


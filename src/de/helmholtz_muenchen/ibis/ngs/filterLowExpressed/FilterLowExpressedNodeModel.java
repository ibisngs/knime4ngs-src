package de.helmholtz_muenchen.ibis.ngs.filterLowExpressed;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;

/**
 * This is the model implementation of FilterLowExpressed.
 * 
 *
 * @author Michael Kluge
 */
public class FilterLowExpressedNodeModel extends RNodeModel {

    // keys for SettingsModels
    protected static final String CFGKEY_KEEP_READS 	= "keepReads";
    protected static final String CFGKEY_KEEP_FRACTION 	= "keepFraction";
    protected static final String CFGKEY_BOTH_SEPERATE	= "bothConditionsSeperate";
    
    // initial default values for SettingsModels
    protected static final int DEFAULT_KEEP_READS = 10;
    protected static final double DEFAULT_KEEP_FRATION = 0.5;
    protected static final boolean DEFAULT_BOTH_SEP = true;
    
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelInteger SET_KEEP_READS	= new SettingsModelInteger(CFGKEY_KEEP_READS, DEFAULT_KEEP_READS);
    private final SettingsModelDouble SET_KEEP_FRACTION	= new SettingsModelDouble(CFGKEY_KEEP_FRACTION, DEFAULT_KEEP_FRATION);
    private final SettingsModelBoolean SET_BOTH_SEP		= new SettingsModelBoolean(CFGKEY_BOTH_SEPERATE, DEFAULT_BOTH_SEP);

    // the logger instance
	@SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FilterLowExpressedNodeModel.class);
	
	private static final String SCRIPT_PATH = "ngs" + File.separatorChar + "de" + File.separatorChar + "filterLowExpressed.R";
	
    /**
     * Constructor for the node model.
     */
	protected FilterLowExpressedNodeModel() {
		super(2, 1, SCRIPT_PATH, new String[]{"--countTable", "--annotationFile"}, new String[]{"--output"});
	}
	
	@Override
	public void init() {
		super.init();
		this.addSetting(SET_KEEP_READS);
		this.addSetting(SET_KEEP_FRACTION);
		this.addSetting(SET_BOTH_SEP);
	}
	
    /**
     * {@inheritDoc}
     * @throws Exception 
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{    	
		BufferedDataTable[] out = super.execute(inData, exec);
		out[0] = exec.createSpecReplacerTable(out[0], this.getSpec(inData[0].getDataTableSpec())); // parse cell types

		return(out);
	}
    
    /**
     * get specs of table
     * @param inSpec
     * @return
     */
    private DataTableSpec getSpec(DataTableSpec inSpec) {
    	return inSpec;
	}
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	
    	// add values set in GUI
    	this.addArgument("--keepReads", this.SET_KEEP_READS.getIntValue());
    	this.addArgument("--keepFraction", this.SET_KEEP_FRACTION.getDoubleValue());
    	this.addArgument("--bothConditionsSeperate", this.SET_BOTH_SEP.getBooleanValue() ? "1" : "0");

    	return new DataTableSpec[]{getSpec(inSpecs[0])};
    }

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void reset() {
		super.reset();
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
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	super.loadInternals(internDir, exec); // load output from stdout and stderr
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	super.saveInternals(internDir, exec); // save output from stdout and stderr
    }
}

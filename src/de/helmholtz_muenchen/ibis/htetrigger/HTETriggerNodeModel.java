package de.helmholtz_muenchen.ibis.htetrigger;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;

/**
 * This is the model implementation of HTETrigger.
 * 
 *
 * @author Tim Jeske
 */
public class HTETriggerNodeModel extends NodeModel {
    
	static final String CFGKEY_USE_HTE = "use_hte";
	static final boolean USE_HTE = true;
	static final String CFGKEY_DEFAULT_THRESHOLD = "threshold";
	static final int DEFAULT_THRESHOLD = 1;
	static final String CFGKEY_LOCAL_THRESHOLD = "local_threshold";
	static final boolean USE_LOCAL_THRESHOLDS = false;
	
	private final SettingsModelBoolean use_hte = new SettingsModelBoolean(HTETriggerNodeModel.CFGKEY_USE_HTE,HTETriggerNodeModel.USE_HTE);
	private final SettingsModelInteger threshold = new SettingsModelInteger(HTETriggerNodeModel.CFGKEY_DEFAULT_THRESHOLD,HTETriggerNodeModel.DEFAULT_THRESHOLD);
	private final SettingsModelBoolean local_threshold = new SettingsModelBoolean(HTETriggerNodeModel.CFGKEY_LOCAL_THRESHOLD,HTETriggerNodeModel.USE_LOCAL_THRESHOLDS);

	int threshold_value, use_hte_value, local_threshold_value;
	
    /**
     * Constructor for the node model.
     */
    protected HTETriggerNodeModel() {
        super(0, 0);
    }

    /**
     * {@inheritDoc}
     */
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
			final ExecutionContext exec){ 

		use_hte_value = 0;
		local_threshold_value = 0;
		threshold_value = DEFAULT_THRESHOLD;
		if (use_hte.getBooleanValue()) {
			use_hte_value = 1;
			threshold_value = threshold.getIntValue();
			if(local_threshold.getBooleanValue()) local_threshold_value = 1;
		}
		pushFlowVariableInt("threshold", threshold_value);
		pushFlowVariableInt("use_hte", use_hte_value);
		pushFlowVariableInt("local_threshold", local_threshold_value);
		return null;
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
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

        return null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         use_hte.saveSettingsTo(settings);
         threshold.saveSettingsTo(settings);
         local_threshold.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	use_hte.loadSettingsFrom(settings);
        threshold.loadSettingsFrom(settings);
        local_threshold.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        use_hte.validateSettings(settings);
        threshold.validateSettings(settings);
        local_threshold.validateSettings(settings);
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


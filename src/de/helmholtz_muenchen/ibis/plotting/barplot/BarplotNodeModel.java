package de.helmholtz_muenchen.ibis.plotting.barplot;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.port.PortObject;
import org.knime.core.node.port.image.ImagePortObject;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RPlottingNodeModel;

/**
 * This is the model implementation of Boxplot.
 * 
 *
 * @author Jonas Zierer
 */
public class BarplotNodeModel extends RPlottingNodeModel {
	public static final String[] PLOT_POINTS_POSSIBILITIES = new String[]{"outliers", "no", "all", "all jittered"};
	
	// the logger instance
    @SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(BarplotNodeModel.class);
      
    // cfgkeys

    // models

    
    /**
     * Constructor for the node model.
     */
    protected BarplotNodeModel() {
    	super(1, "plotting" + File.separatorChar + "barplot.R", new String[]{"--data"}, new String[]{"--output"});
    }

    /**
     * {@inheritDoc}
     * @throws CanceledExecutionException 
     */
    @Override
    protected ImagePortObject[] execute(final PortObject[] inObjects, final ExecutionContext exec) throws Exception {
    	
        return(super.execute(inObjects, exec));
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
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	super.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	super.loadInternals(internDir, exec);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	super.saveInternals(internDir, exec);
    }
}


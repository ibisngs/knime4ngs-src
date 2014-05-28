package de.helmholtz_muenchen.ibis.plotting.histogram;

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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.port.PortObject;
import org.knime.core.node.port.image.ImagePortObject;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RPlottingNodeModel;

/**
 * This is the model implementation of Boxplot.
 * 
 *
 * @author Jonas Zierer
 */
public class HistogramNodeModel extends RPlottingNodeModel {
	public static final String[] PLOT_POINTS_POSSIBILITIES = new String[]{"outliers", "no", "all", "all jittered"};
	
	// the logger instance
    @SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(HistogramNodeModel.class);
      
    // cfgkeys
    static public final String CFGKEY_DENSITY         = "density";
    static public final String CFGKEY_DENSITY_CURVE   = "density curve";
    static public final String CFGKEY_BINWIDTH        = "binwidth";
	
    // models
    protected final SettingsModelBoolean m_density       = new SettingsModelBoolean(HistogramNodeModel.CFGKEY_DENSITY, false);
    protected final SettingsModelBoolean m_density_curve = new SettingsModelBoolean(HistogramNodeModel.CFGKEY_DENSITY_CURVE, false);
    protected final SettingsModelDouble m_binwidth       = new SettingsModelDouble(HistogramNodeModel.CFGKEY_BINWIDTH, 1.0);
    
    
    /**
     * Constructor for the node model.
     */
    protected HistogramNodeModel() {
    	super(1, "plotting" + File.separatorChar + "histogram.R", new String[]{"--data"}, new String[]{"--output"});
    }

    /**
     * {@inheritDoc}
     * @throws CanceledExecutionException 
     */
    @Override
    protected ImagePortObject[] execute(final PortObject[] inObjects, final ExecutionContext exec) throws CanceledExecutionException {
    	this.setFlag("--dens"      , m_density.getBooleanValue());
    	this.setFlag("--densCurve" , m_density_curve.getBooleanValue());
    	this.addArgument("--binwidth", m_binwidth.getDoubleValue());
    	
    	this.m_col_y.setEnabled(false);
    	
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
    	m_density.saveSettingsTo(settings);
    	m_density_curve.saveSettingsTo(settings);
    	m_binwidth.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
    	m_density.loadSettingsFrom(settings);
    	m_density_curve.loadSettingsFrom(settings);
    	m_binwidth.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.validateSettings(settings);
    	m_density.validateSettings(settings);
    	m_density_curve.validateSettings(settings);
    	m_binwidth.validateSettings(settings);
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


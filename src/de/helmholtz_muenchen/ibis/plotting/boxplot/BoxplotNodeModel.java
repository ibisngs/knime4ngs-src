package de.helmholtz_muenchen.ibis.plotting.boxplot;

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
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.port.PortObject;
import org.knime.core.node.port.image.ImagePortObject;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RPlottingNodeModel;

/**
 * This is the model implementation of Boxplot.
 * 
 *
 * @author Jonas Zierer
 */
public class BoxplotNodeModel extends RPlottingNodeModel {
	public static final String[] PLOT_POINTS_POSSIBILITIES = new String[]{"outliers", "no", "all", "all jittered"};
	
	// the logger instance
    @SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(BoxplotNodeModel.class);
      
    // cfgkeys
    static public final String CFGKEY_POINTS             = "plot points";
    static public final String CFGKEY_BOXPERCOL          = "one box per column";
    static public final String CFGKEY_BOXPERCOL_COLUMNS  = "columns selection";
	
    // models
    protected final SettingsModelString m_points                  = new SettingsModelString(BoxplotNodeModel.CFGKEY_POINTS, BoxplotNodeModel.PLOT_POINTS_POSSIBILITIES[0]);
    protected final SettingsModelBoolean m_boxpercol              = new SettingsModelBoolean(BoxplotNodeModel.CFGKEY_BOXPERCOL, false);
    protected final SettingsModelFilterString m_boxpercol_columns = new SettingsModelFilterString(BoxplotNodeModel.CFGKEY_BOXPERCOL_COLUMNS);
   
    /**
     * Constructor for the node model.
     */
    protected BoxplotNodeModel() {
    	super(1, "plotting" + File.separatorChar + "boxplot.R", new String[]{"--data"}, new String[]{"--output"});
    }

    /**
     * {@inheritDoc}
     * @throws CanceledExecutionException 
     */
    @Override
    protected ImagePortObject[] execute(final PortObject[] inObjects, final ExecutionContext exec) throws CanceledExecutionException {
    	this.addArgument("--points", m_points.getStringValue());
    	if(m_boxpercol.getBooleanValue()){
    		this.addArgument("--columns", m_boxpercol_columns.getIncludeList());
    		this.m_col_x.setEnabled(false);
    		this.m_col_y.setEnabled(false);
    	}else{
    		this.m_col_x.setEnabled(true);
    		this.m_col_y.setEnabled(true);
    	}
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
    	m_points.saveSettingsTo(settings);
    	m_boxpercol.saveSettingsTo(settings);
    	m_boxpercol_columns.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
    	m_points.loadSettingsFrom(settings);
    	m_boxpercol.loadSettingsFrom(settings);
    	m_boxpercol_columns.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.validateSettings(settings);
    	m_points.validateSettings(settings);
    	m_boxpercol.validateSettings(settings);
    	m_boxpercol_columns.validateSettings(settings);
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


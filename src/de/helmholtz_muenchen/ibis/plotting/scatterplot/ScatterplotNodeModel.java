package de.helmholtz_muenchen.ibis.plotting.scatterplot;

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
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.port.PortObject;
import org.knime.core.node.port.image.ImagePortObject;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RPlottingNodeModel;

/**
 * This is the model implementation of Boxplot.
 * 
 *
 * @author Jonas Zierer
 */
public class ScatterplotNodeModel extends RPlottingNodeModel {
	// the logger instance
    @SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(ScatterplotNodeModel.class);
      
    // cfgkeys
    static public final String CFGKEY_MATRIX          = "scattermatrix";
    static public final String CFGKEY_MATRIX_COLUMNS  = "scattermatrix columns";
    static public final String CFGKEY_ALPHA           = "alpha";
    static public final String CFGKEY_POINTSIZE       = "pointsize";
    
    // models
    protected final SettingsModelBoolean m_matrix = new SettingsModelBoolean(ScatterplotNodeModel.CFGKEY_MATRIX, false);
    private final SettingsModelFilterString m_matrix_columns   = new SettingsModelFilterString(ScatterplotNodeModel.CFGKEY_MATRIX_COLUMNS);
    protected final SettingsModelDoubleBounded m_alpha = new SettingsModelDoubleBounded(ScatterplotNodeModel.CFGKEY_ALPHA, 0.8, 0.0, 1.0);
    protected final SettingsModelIntegerBounded m_pointsize = new SettingsModelIntegerBounded(ScatterplotNodeModel.CFGKEY_POINTSIZE, 1, 1, 20);

    
    /**
     * Constructor for the node model.
     */
    protected ScatterplotNodeModel() {
    	super(1, "plotting" + File.separatorChar + "scatterplot.R", new String[]{"--data"}, new String[]{"--output"});
    }

    /**
     * {@inheritDoc}
     * @throws CanceledExecutionException 
     */
    @Override
    protected ImagePortObject[] execute(final PortObject[] inObjects, final ExecutionContext exec) throws CanceledExecutionException {
    	if(m_matrix.getBooleanValue()){
    		this.addArgument("--matrix", m_matrix_columns.getIncludeList());
    		this.m_col_x.setEnabled(false);
    		this.m_col_y.setEnabled(false);
    	}else{
    		this.m_col_x.setEnabled(true);
    		this.m_col_y.setEnabled(true);
    	}
    	this.addArgument("--alpha", m_alpha.getDoubleValue());
    	this.addArgument("--size", m_pointsize.getIntValue());
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
    	m_matrix.saveSettingsTo(settings);
    	m_matrix_columns.saveSettingsTo(settings);
    	m_alpha.saveSettingsTo(settings);
    	m_pointsize.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
    	m_matrix.loadSettingsFrom(settings);
    	m_matrix_columns.loadSettingsFrom(settings);
    	m_alpha.loadSettingsFrom(settings);
    	m_pointsize.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.validateSettings(settings);
    	m_matrix.validateSettings(settings);
    	m_matrix_columns.validateSettings(settings);
    	m_alpha.validateSettings(settings);
    	m_pointsize.validateSettings(settings);
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


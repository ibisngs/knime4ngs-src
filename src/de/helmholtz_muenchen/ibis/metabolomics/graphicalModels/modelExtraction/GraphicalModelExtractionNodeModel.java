package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.modelExtraction;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * @author Jonas Zierer
 */
public class GraphicalModelExtractionNodeModel extends RNodeModel {

	/** CFG KEYS */
	static final String CFGKEY_PERC_INCLUSION = "percentageInclusion";
	static final String CFGKEY_EV    = "e.v.";

    
	/** SETTING MODELS */
	private final SettingsModelDoubleBounded  m_percIncl = new SettingsModelDoubleBounded(GraphicalModelExtractionNodeModel.CFGKEY_PERC_INCLUSION, 0.8, 0.0, 1.0);
    private final SettingsModelIntegerBounded m_ev       = new SettingsModelIntegerBounded(GraphicalModelExtractionNodeModel.CFGKEY_EV, 5, 0, Integer.MAX_VALUE);

    
    /**
     * Constructor for the node model.
     */
    protected GraphicalModelExtractionNodeModel() {
        super(1, 2, "statistics" + File.separatorChar + "graphicalModels" + File.separatorChar + "runner.modelExtraction.R", new String[]{"--input"}, new String[]{"--edgeslist", "--adjacency"});
    }

    /**
     * {@inheritDoc}
     * @throws Exception 
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{
    	this.addArgument("--percIncl", m_percIncl.getDoubleValue());
    	this.addArgument("--ev"      , m_ev.getIntValue());

    	
		return(super.execute(inData, exec));
	}
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
        return new DataTableSpec[]{null, null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
        m_ev.saveSettingsTo(settings);
        m_percIncl.saveSettingsTo(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
        m_ev.loadSettingsFrom(settings);
        m_percIncl.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	m_ev.validateSettings(settings);
        m_percIncl.validateSettings(settings);

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


package de.helmholtz_muenchen.ibis.misc.sleepnode;

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
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;


/**
 * This is the model implementation of SleepNode.
 * Node sleeps for the given time
 *
 * @author Jonas Zierer
 */
public class SleepNodeNodeModel extends NodeModel {
    
	public static final String[] POSSIBLE_STRING_VALUES = {"String 1", "String 2", "Another String"};
	
	/** CFG KEYS */
	static final String CFGKEY_TIME    = "argInt";
    
	/** DEFAULTS */
	static final int DEFAULT_TIME = 100;
	/** SETTING MODELS */
    private final SettingsModelIntegerBounded m_time  = new SettingsModelIntegerBounded(SleepNodeNodeModel.CFGKEY_TIME, 5, 0, 1000000);

    
    /**
     * Constructor for the node model.
     */
    protected SleepNodeNodeModel() {
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     * @throws Exception 
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{
		Thread.sleep(this.m_time.getIntValue());
	
		return(inData);
	}
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
        return inSpecs;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
        m_time.saveSettingsTo(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
        m_time.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	m_time.validateSettings(settings);

    }

	@Override
	protected void loadInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {
		
	}

	@Override
	protected void saveInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {		
	}

	@Override
	protected void reset() {
	}

}


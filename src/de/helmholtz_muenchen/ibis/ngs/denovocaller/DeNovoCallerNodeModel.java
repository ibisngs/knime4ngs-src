package de.helmholtz_muenchen.ibis.ngs.denovocaller;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


/**
 * This is the model implementation of DeNovoCaller.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class DeNovoCallerNodeModel extends NodeModel {
    
	/**
	 * Config Keys
	 */
	public static final String CFGKEY_INPUTTYPE = "inputtype";
	public static final String CFGKEY_PEDFILE = "pedfile";
	
	
	
	/**
	 * Node Models
	 */
	private final SettingsModelString m_inputtype = new SettingsModelString(DeNovoCallerNodeModel.CFGKEY_INPUTTYPE,"");
	private final SettingsModelString m_pedfile = new SettingsModelString(DeNovoCallerNodeModel.CFGKEY_PEDFILE,"");
	
	/**
	 * Logger
	 */
	private static final NodeLogger LOGGER = NodeLogger.getLogger(DeNovoCallerNodeModel.class);
	
    /**
     * Constructor for the node model.
     */
    protected DeNovoCallerNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String inputtype = "";
    	if(m_inputtype.getStringValue().equals("Multi-Sample VCFs (Trios)")){
    		LOGGER.info("Using 'Trio' filtering.");
    		inputtype = "Trio";
    	}else{
    		LOGGER.info("Using 'Single' filtering.");
    		inputtype = "Single";
    	}
    	
    	BufferedDataTable[] outTable = DeNovoCaller.findDeNovos(inData[0], m_pedfile.getStringValue(),inputtype, exec, LOGGER);
    	
        return outTable;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

        // TODO: generated method stub
        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
        m_inputtype.saveSettingsTo(settings);
        m_pedfile.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
         m_inputtype.loadSettingsFrom(settings);
         m_pedfile.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_inputtype.validateSettings(settings);
        m_pedfile.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }

}


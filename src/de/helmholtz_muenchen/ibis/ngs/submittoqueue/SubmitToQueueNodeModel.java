package de.helmholtz_muenchen.ibis.ngs.submittoqueue;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

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
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * This is the model implementation of SubmitToQueue.
 * 
 *
 * @author 
 */
public class SubmitToQueueNodeModel extends NodeModel {
	
	public static final String CFGKEY_MEMORY = "memory";
	public static final String CFGKEY_JOBNAME = "jobname";

	public static final int DEFAULT_MEMORY = 5;

	private final SettingsModelIntegerBounded m_memory = new SettingsModelIntegerBounded(CFGKEY_MEMORY,DEFAULT_MEMORY,1,100);
	private final SettingsModelString m_jobname = new SettingsModelString(CFGKEY_JOBNAME, "");
	
	
    /**
     * Constructor for the node model.
     */
    protected SubmitToQueueNodeModel() {
    	
        super(0, 0);
        
        m_jobname.setStringValue("NGS");
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	int memory = m_memory.getIntValue();
    	System.out.println(memory);
    	String jobPreName = m_jobname.getStringValue(); 
    	
    	pushFlowVariableString("Memory", memory+"");
    	pushFlowVariableString("JobPrefix", jobPreName);

        return new BufferedDataTable[]{};
        
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
    	
    	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", "qstat -h");
    	Process p = null;
		try {
			p = b.start();
		} catch (IOException e) {
			e.printStackTrace();
		}
    	try {
			p.waitFor();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
    	BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getErrorStream()));
    	String s = null;
    	StringBuffer nodeEntry = new StringBuffer(60);
    	try {
			while ((s = stdInput.readLine()) != null) {
				nodeEntry.append(s+" ");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
    	int pos = nodeEntry.toString().lastIndexOf("job");
    	if(pos == -1) {
    		throw new InvalidSettingsException("There is no Sun Grid Engine installed on your system. Unfortunately you cannot use this node. Possibly you can try to install or configure a queue on your system.");
    	}
    	
        return new DataTableSpec[]{};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_jobname.saveSettingsTo(settings);
    	m_memory.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_jobname.loadSettingsFrom(settings);
    	m_memory.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_jobname.validateSettings(settings);
    	m_memory.validateSettings(settings);
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


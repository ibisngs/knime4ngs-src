package de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode;

import java.io.File;
import java.io.IOException;

import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;

import de.helmholtz_muenchen.ibis.utils.Utils;

public abstract class ScriptNodeModel extends NodeModel {
	protected final StringBuffer STDOUT;
	protected final StringBuffer STDERR;
	protected final String SCRIPT;

	private static final NodeLogger LOGGER = NodeLogger.getLogger(ScriptNodeModel.class);


	protected ScriptNodeModel(int nrInDataPorts, int nrOutDataPorts, String script) {
		super(nrInDataPorts, nrOutDataPorts);
		this.STDERR = new StringBuffer();
		this.STDOUT = new StringBuffer();
		
		this.SCRIPT = getScriptPath() + script;
	}
	
	protected String getScriptPath(){
		LOGGER.error("CALLING SCRIPT PATH METHOD OF SCRIPT NODE MODEL");
		return("");
	}
	
	protected void executeScript(String args, final ExecutionContext exec) throws CanceledExecutionException {
		this.STDERR.setLength(0);
		this.STDOUT.setLength(0);
		
		// run command
		String com= this.SCRIPT + args;

		LOGGER.info("RUNNING SCRIPT: >" + com);
		try {
			Utils.executeCommand(com, this.STDOUT, this.STDERR, LOGGER);
		} catch (Exception e) {
			//LOGGER.error("Couldn'r run command >" + com);
			//LOGGER.error("STDERR: " + this.STDERR + "\n\nSTDOUT: " + this.STDOUT);
			throw new CanceledExecutionException("Error occured while ececuting external command >" + com + "\n" + e.getMessage());
		}
		
	}

	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// OVERIDE KNIME NODE METHODS
	/////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
    	this.STDERR.setLength(0);
		this.STDOUT.setLength(0);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
        // TODO load STDOUT and STDERR
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	// TODO save STDOUT and STDERR
    }

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// GETTERS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	public StringBuffer getSTDOUT() {
		return STDOUT;
	}
	public StringBuffer getSTDERR() {
		return STDERR;
	}
	public String getSCRIPT() {
		return SCRIPT;
	}
	

}

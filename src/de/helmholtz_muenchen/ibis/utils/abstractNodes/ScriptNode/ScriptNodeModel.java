package de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.ExecutionException;

import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.port.PortType;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode.ExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;

public abstract class ScriptNodeModel extends ExecutorNodeModel {
	protected final String SCRIPT;
	public static final String SCRIPTS_SUBFOLDER = "scripts";
	
	
	protected ScriptNodeModel(int nrInDataPorts, int nrOutDataPorts, String script) {
		super(nrInDataPorts, nrOutDataPorts, true, true);		
		this.SCRIPT = getScriptPath() + script;
	}
	
	protected ScriptNodeModel(final PortType[] inPortTypes, final PortType[] outPortTypes, String script) {
		super(inPortTypes, outPortTypes, true, true);		
		this.SCRIPT = getScriptPath() + script;
	}
	
	protected String getScriptPath(){
		return(IO.getScriptPath() + SCRIPTS_SUBFOLDER + File.separatorChar);
	}
	
	protected void executeScript(final ExecutionContext exec, String[] environment) throws CanceledExecutionException, IOException, InterruptedException, ExecutionException, UnsuccessfulExecutionException {
		executeCommand(exec, this.getCommand(), environment, null, null);
	}
	
	protected abstract String[] getCommand();
	

	public String getSCRIPT() {
		return SCRIPT;
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// OVERIDE KNIME NODE METHODS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void reset() {
		super.reset();
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

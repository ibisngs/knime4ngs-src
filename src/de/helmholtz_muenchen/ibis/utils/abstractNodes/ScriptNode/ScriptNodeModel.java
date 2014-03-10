package de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode;

import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode.ExecutorNodeModel;

public abstract class ScriptNodeModel extends ExecutorNodeModel {
	protected final String SCRIPT;

	protected ScriptNodeModel(int nrInDataPorts, int nrOutDataPorts, String script) {
		super(nrInDataPorts, nrOutDataPorts, true, true);		
		this.SCRIPT = getScriptPath() + script;
	}
	
	protected String getScriptPath(){
		return(IO.getScriptPath());
	}
	
	protected void executeScript(final ExecutionContext exec, String[] environment) throws CanceledExecutionException {
		executeCommand(exec, this.getCommand(), environment, true);
	}
	
	protected abstract String[] getCommand();
	

	public String getSCRIPT() {
		return SCRIPT;
	}
}

package de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

public abstract class ScriptNodeModel extends NodeModel {
	protected String STDOUT;
	protected String STDERR;
	protected final String SCRIPT;

	private static final NodeLogger LOGGER = NodeLogger.getLogger(ScriptNodeModel.class);


	protected ScriptNodeModel(int nrInDataPorts, int nrOutDataPorts, String script) {
		super(nrInDataPorts, nrOutDataPorts);
		this.STDERR = "";
		this.STDOUT = "";
		
		this.SCRIPT = getScriptPath() + script;
	}
	
	protected String getScriptPath(){
		return(IO.getScriptPath());
	}
	
	protected void executeScript(final ExecutionContext exec, String[] environment) throws CanceledExecutionException {
		try {
			StringBuffer stdout = new StringBuffer();
			StringBuffer stderr = new StringBuffer();
			Executor.executeCommand(this.getCommand(), exec, environment, LOGGER, stdout, stderr);
			this.STDERR = stderr.toString();
			this.STDOUT = stdout.toString();
		} catch (Exception e) {
			LOGGER.error(e.getMessage());
			throw(new CanceledExecutionException(e.getMessage()));
		}
	}

	protected abstract String[] getCommand();
	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// OVERIDE KNIME NODE METHODS
	/////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
    	this.STDERR = "";
		this.STDOUT = "";
    }
    
    public static final String FILE_STDOUT = "stdout.txt";
    public static final String FILE_STDERR = "stderr.txt";
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	File f_stdout = new File(internDir, FILE_STDOUT);
    	File f_stderr = new File(internDir, FILE_STDERR);
        
    	if(f_stdout.exists()){
    		this.STDOUT = (FileUtils.readFileToString(f_stdout));
    	}
    	if(f_stderr.exists()){
    		this.STDERR = (FileUtils.readFileToString(f_stderr));
    	}
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	File f_stdout = new File(internDir, FILE_STDOUT);
    	File f_stderr = new File(internDir, FILE_STDERR);
    	
    	FileUtils.write(f_stdout, this.STDOUT);
    	FileUtils.write(f_stderr, this.STDERR);
    }

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// GETTERS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	public String getSTDOUT() {
		return STDOUT;
	}
	public String getSTDERR() {
		return STDERR;
	}
	public String getSCRIPT() {
		return SCRIPT;
	}
	

}

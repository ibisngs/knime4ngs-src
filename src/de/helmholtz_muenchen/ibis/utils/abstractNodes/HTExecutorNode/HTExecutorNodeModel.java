package de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode;

import java.io.File;
import java.io.IOException;
import java.net.InetAddress;
import java.util.Map;

import javax.swing.JOptionPane;

import org.knime.core.node.*;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.port.PortType;
import org.knime.core.node.workflow.FlowVariable;

import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.threads.ExecuteThread;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;

/**
 * This is the model implementation of HTExecutorNode.
 * 
 * 
 * @author Tim Jeske
 */
public abstract class HTExecutorNodeModel extends NodeModel {
	
	//variables characterizing the HTE_Node
	private int hte_id = -1;
	private int node_id = -1;
	private String lockCommand = "";
	private String host_name = "";
	private String node_name = this.getClass().getCanonicalName();
	private String defaultLockFile = "/home/ibis/tim.jeske/NGSTest/"+node_name+SuccessfulRunChecker.LOCK_ENDING;
	private int count = 0;
	private int threshold_value = DEFAULT_THRESHOLD;
	
	//internal flags
	private boolean use_hte_value = false;
	private int resetCounter = 0;
	private boolean calledByException = false;
	private boolean executionTried = false;
	
	static final String CFGKEY_DEFAULT_THRESHOLD = "threshold";
	static final int DEFAULT_THRESHOLD = 1;

	private final SettingsModelInteger threshold = new SettingsModelInteger(
			HTExecutorNodeModel.CFGKEY_DEFAULT_THRESHOLD, 1);

	protected HTExecutorNodeModel(PortType[] inPortTypes,
			PortType[] outPortTypes) {
		super(inPortTypes, outPortTypes);
	}

	protected HTExecutorNodeModel(int nrInDataPorts, int nrOutDataPorts) {
		super(nrInDataPorts, nrOutDataPorts);
	}

	private void readFlowVariables(){
		
		Map<String, FlowVariable> map = this.getAvailableFlowVariables();

		if (!map.containsKey("use_hte")) {
			System.err
			.println("No information about HTE usage found! Please add the HTETrigger, if you want to use HTE.");
			use_hte_value=false;
			return;
		}
		
		if (map.get("use_hte").getIntValue() == 1) {
			use_hte_value = true;
		} else {
			use_hte_value = false;
		}

		if (map.containsKey("hte_id")) {
			hte_id = map.get("hte_id").getIntValue();
		} 
		
		boolean use_local_threshold = false;
		
		if (map.containsKey("local_threshold")) {
			use_local_threshold = (map.get("local_threshold").getIntValue()==1);
		}
		
		if(use_local_threshold) {
			threshold_value = threshold.getIntValue();
		} else {
			threshold_value = map.get("threshold").getIntValue();
		}
	}
	
	private void recExecuteCommand(String[] command, ExecutionContext exec,
			String[] environment, NodeLogger logger, String stdOutFile,
			String stdErrFile, StringBuffer stdOut, StringBuffer stdErr,
			String StdInFile, SuccessfulRunChecker checker, HTEDBHandler htedb) throws Exception {

		int exitcode = 0;
		
		if (count < threshold_value) {
			count++;
			htedb.updateCount(hte_id, node_id, count);
			exitcode = Executor.executeCommandWithExitCode(command, exec,
					environment, logger, stdOutFile, stdErrFile, stdOut,
					stdErr, StdInFile);
			if (exitcode == 0) {
				htedb.writeSuccess(hte_id, node_id);
				checker.writeOK();
				checker.finalize();
				return;
			} else if (exitcode == 139) {
				htedb.writeError(hte_id, node_id, "segmentation fault");
			} else {
				htedb.writeError(hte_id, node_id, "exitcode: " + exitcode);
			}
		}

		if (count >= threshold_value) {
			System.err.println("Email: Your command " + command[0]
					+ " could not be completed successfully!");
			
			htedb.closeConnection();
			calledByException = true;
			throw (new UnsuccessfulExecutionException("Exit code was not 0: '"
					+ exitcode + "' for " + ExecuteThread.getCommand(command)));

		} else {
			recExecuteCommand(command, exec, environment, logger, stdOutFile,
					stdErrFile, stdOut, stdErr, StdInFile, checker, htedb);
		}
	}

	protected void executeCommand(String[] command, ExecutionContext exec,
			String[] environment, NodeLogger logger,File lockFile, String stdOutFile,
			String stdErrFile, StringBuffer stdOut, StringBuffer stdErr,
			String StdInFile) throws Exception {

		calledByException = false;
		
		this.readFlowVariables();

		//prepare klock
		host_name = InetAddress.getLocalHost().getHostName();
		lockCommand = "";
		for (String s : command) {
			lockCommand += s;
		}

		logger.info(node_name+" is executed with: use_hte="+use_hte_value+" hte_id="+hte_id+" threshold="+threshold_value);
		
		if(lockFile == null) {
			lockFile = new File(defaultLockFile);
		}
		boolean terminationState = SuccessfulRunChecker
				.hasTerminatedSuccessfully(lockFile, lockCommand);
		logger.info("Successful termination state: " + terminationState);

		//abort execution if node has been executed successfully
		if (terminationState) {
			System.out.println("Node "+node_name+" has been executed succesfully according to klock!");
			return;
		}
			
		SuccessfulRunChecker checker = new SuccessfulRunChecker(lockFile,
				lockCommand);
		
		//if HTE is used
		if (use_hte_value) {
			if(hte_id ==-1) {
				throw new UnsuccessfulExecutionException("Execute HTETrigger Node before the other nodes in the workflow!");
			}
			HTEDBHandler htedb = new HTEDBHandler();
			node_id = htedb.getHTENodeId(lockCommand, node_name, host_name, hte_id);
			System.out.println("Threshold: "+threshold_value);
			if(node_id == -1) {
				node_id = htedb.insertNewHTENode(lockCommand, node_name, host_name,threshold_value, hte_id);
			} else {
				int tmp_threshold = htedb.getHTENodeThreshold(node_id, hte_id);
				if(tmp_threshold!=threshold_value) {
					htedb.updateThreshold(hte_id, node_id, threshold_value);
				}
			}
			if(count>=threshold_value) {
				this.popUpWindow(htedb);
			}
			executionTried = true;
			recExecuteCommand(command, exec, environment, logger, stdOutFile,
					stdErrFile, stdOut, stdErr, StdInFile, checker, htedb);
			htedb.closeConnection();
		} else {
			//HTE is not used
			logger.info("HTE is not used");
			Executor.executeCommand(command, exec, environment, logger,
					stdOutFile, stdErrFile, stdOut, stdErr, StdInFile);
			checker.writeOK();
			checker.finalize();
		}
	}
	
	private void popUpWindow(HTEDBHandler htedb) {
		String newline = System.getProperty("line.separator");
		String message = node_name + " was executed " + count + " times."
				+ newline + "Your threshold was " + threshold_value+ "." + newline
				+ "Do you want to reset the count in the DB?";
		int answer = JOptionPane.showOptionDialog(null, message, "Reset in database?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, null, null);
		if(answer == JOptionPane.YES_OPTION) {
			count = 0;
			
			htedb.updateCount(hte_id, node_id, 0);
		}	
	}

	/**
	 * called before reset operations of knime node
	 */
	@Override
	protected void reset() {
		if(!executionTried) return;
		
		if(calledByException) {
			resetCounter++;
			if(resetCounter<=2) return;
			calledByException = false;
			resetCounter = 0;
		}
		
		this.readFlowVariables();
		
		if(use_hte_value) {
			HTEDBHandler htedb = new HTEDBHandler();
			this.popUpWindow(htedb);
			htedb.closeConnection();
		}
		
		/*
		 
	 	//should only be executed if reset is called explicitly 
		
		Map<String, FlowVariable> map = this.getAvailableFlowVariables();
		
		//return if HTE is off
		if (!map.containsKey("use_hte")) return;
		if (map.get("use_hte").getIntValue() == 0) return;
		
		//HTE is on
		htedb = new HTEDBHandler();
		this.popUpWindow();
		
		if (map.containsKey("hte_id")) {
			hte_id = map.get("hte_id").getIntValue();
		}
		
		node_id = htedb.getHTENodeId(lockCommand, node_name, host_name, hte_id);
		if(node_id == -1) return;
		
		
    	threshold_value = htedb.getHTENodeThreshold(node_id,hte_id);
    	
		count = htedb.getHTENodeCount(node_id,hte_id);
		
		//HTENode entry is in database
		this.popUpWindow();
		htedb.closeConnection();
		
		*/
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveSettingsTo(final NodeSettingsWO settings) {
		threshold.saveSettingsTo(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
			throws InvalidSettingsException {
		threshold.loadSettingsFrom(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void validateSettings(final NodeSettingsRO settings)
			throws InvalidSettingsException {
		threshold.validateSettings(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadInternals(final File internDir,
			final ExecutionMonitor exec) throws IOException,
			CanceledExecutionException {
    	/*File f_host = new File(internDir, "host.txt");
    	File f_lockCommand = new File(internDir, "command.txt");
        
    	// set the buffers, if the files are there
    	if(f_host.exists() && f_lockCommand.exists()){

    		host_name = FileUtils.readFileToString(f_host);
    		lockCommand= FileUtils.readFileToString(f_lockCommand);

    	}*/
    	
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveInternals(final File internDir,
			final ExecutionMonitor exec) throws IOException,
			CanceledExecutionException {
		/*File f_host = new File(internDir, "host.txt");
    	File f_lockCommand = new File(internDir, "command.txt");
    	
    	FileUtils.write(f_host, host_name);
    	FileUtils.write(f_lockCommand, lockCommand);*/
	}

}

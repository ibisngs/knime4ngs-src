package de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode;

import java.io.File;
import java.net.InetAddress;
import java.sql.SQLException;

import org.knime.core.node.*;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.port.PortType;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
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
	
	//variables characterizing the HTExecution
	private int exec_id = -1;
	private String lockCommand = "";
	private String host_name = "";
	private String node_name = this.getClass().getCanonicalName();
	private String defaultLockFile = System.getProperty("user.home")+"/"+node_name+SuccessfulRunChecker.LOCK_ENDING;
	private int count = 0;
	private int threshold_value = DEFAULT_THRESHOLD;
	private String db_file = "";
	
	//internal flags
	private boolean use_hte = false;
	
	static final String CFGKEY_DEFAULT_THRESHOLD = "threshold";
	static final int DEFAULT_THRESHOLD = 1;
	
	static final String CFGKEY_USE_PREF = "use_pref";

	
	private final SettingsModelInteger threshold = new SettingsModelInteger(
			HTExecutorNodeModel.CFGKEY_DEFAULT_THRESHOLD, 1);

	protected HTExecutorNodeModel(PortType[] inPortTypes,
			PortType[] outPortTypes) {
		super(inPortTypes, outPortTypes);
	}

	protected HTExecutorNodeModel(int nrInDataPorts, int nrOutDataPorts) {
		super(nrInDataPorts, nrOutDataPorts);
	}
	
	private void recExecuteCommand(String[] command, ExecutionContext exec,
			String[] environment, NodeLogger logger, String stdOutFile,
			String stdErrFile, StringBuffer stdOut, StringBuffer stdErr,
			String StdInFile, SuccessfulRunChecker checker, HTEDBHandler htedb) throws Exception {

		int exitcode = 0;
		
		if (count < threshold_value) {
			count++;
			htedb.updateCount(exec_id, count);
			exitcode = Executor.executeCommandWithExitCode(command, exec,
					environment, logger, stdOutFile, stdErrFile, stdOut,
					stdErr, StdInFile);
			if (exitcode == 0) {
				htedb.writeSuccess(exec_id);
				checker.writeOK();
				checker.finalize();
				return;
			} else if (exitcode == 139) {
				htedb.writeError(exec_id, "segmentation fault");
			} else {
				htedb.writeError(exec_id, "exitcode: " + exitcode);
			}
		}

		if (count >= threshold_value) {
			System.err.println("Email: Your command " + command[0]
					+ " could not be completed successfully!");
			
			htedb.closeConnection();
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
		
		use_hte = IBISKNIMENodesPlugin.getDefault().getHTEPreference();
		threshold_value = threshold.getIntValue();
		db_file = IBISKNIMENodesPlugin.getDefault().getDBFilePreference();

		this.count = 0;
		//prepare klock
		host_name = InetAddress.getLocalHost().getHostName();
		lockCommand = "";
		for (String s : command) {
			lockCommand += s;
		}

		logger.info(node_name+" is executed with: use_hte="+use_hte+" threshold="+threshold_value);
		
		if(lockFile == null) {
			lockFile = new File(defaultLockFile);
		}
		logger.info("klock file can be found in "+lockFile);
		boolean terminationState = SuccessfulRunChecker
				.hasTerminatedSuccessfully(lockFile, lockCommand);
		logger.info("Successful termination state: " + terminationState);

		//abort execution if node has been executed successfully
		if (terminationState) {
			return;
		}
			
		SuccessfulRunChecker checker = new SuccessfulRunChecker(lockFile,
				lockCommand);
		
		HTEDBHandler htedb = null;
		
		//try to establish connection to database
		if (use_hte) {
			try {
				htedb = new HTEDBHandler(db_file);
			} catch (SQLException e) {
				System.err.println("Connection to database could not be established: "+e.getMessage());
				use_hte = false;
			}
		}
		
		if (use_hte) {
			exec_id = htedb.insertNewHTExecution(lockCommand, node_name, host_name, threshold_value);
//			System.out.println("########################## this is the exec_id: "+exec_id);
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
}

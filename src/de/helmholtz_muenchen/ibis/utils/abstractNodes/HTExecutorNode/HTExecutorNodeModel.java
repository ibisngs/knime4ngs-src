package de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode;

import java.io.File;
import java.net.InetAddress;
import java.sql.SQLException;
import java.util.Properties;

import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;

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
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(HTExecutorNodeModel.class);
	
	//variables characterizing the HTExecution
	private int exec_id = -1;
	private String lockCommand = "";
	private String host_name = "";
	private String node_name = this.getClass().getCanonicalName();
	private String defaultLockFile = System.getProperty("user.home")+"/"+node_name+SuccessfulRunChecker.LOCK_ENDING;
	private int count = 0;
	private int threshold_value = DEFAULT_THRESHOLD;
	private String db_file = "";
	private String email = "";
	
	private final String EMAILSENDER = "ibis.knime@helmholtz-muenchen.de";
	private final String EMAILHOST = "outmail.helmholtz-muenchen.de";
	private final String HEADER = "IBIS KNIME Nodes Notification";
	
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
			String[] environment, String stdOutFile,
			String stdErrFile, StringBuffer stdOut, StringBuffer stdErr,
			String StdInFile, SuccessfulRunChecker checker, HTEDBHandler htedb) throws Exception {

		int exitcode = 0;
		String err_msg = "";
		
		if (count < threshold_value) {
			count++;
			htedb.updateCount(exec_id, count);
			exitcode = Executor.executeCommandWithExitCode(command, exec,
					environment,LOGGER, stdOutFile, stdErrFile, stdOut,
					stdErr, StdInFile);
			err_msg = stdErr.toString();
			stdErr = new StringBuffer();
			if (exitcode == 0) {
				htedb.writeSuccess(exec_id);
				checker.writeOK();
				checker.finalize();
				return;
			} else {
				htedb.writeError(exec_id, err_msg);
			}
		}

		if (count >= threshold_value) {
			String cmd = "";
			if(command.length==1) {
				cmd = " " +command[0];
			} else {
				for(String c : command) {
					cmd = " "+ c;
				}
			}
			
			String newLine = System.getProperty("line.separator");
			
			if(email!=null) {
				sendMail("Hello " + email.split("@")[0] + ", "
						+ newLine
						+ "unfortunately the node " + node_name + " failed " + count
						+ " time(s) on host " + host_name + "."
						+ newLine + newLine 
						+ "Your command was: "
						+ newLine
						+ cmd 
						+ newLine + newLine
						+ "The error message was: "
						+ newLine
						+ err_msg);
			}
			
			htedb.closeConnection();
			throw (new UnsuccessfulExecutionException("Exit code was not 0: '"
					+ exitcode + "' for " + ExecuteThread.getCommand(command)));

		} else {
			recExecuteCommand(command, exec, environment, stdOutFile,
					stdErrFile, stdOut, stdErr, StdInFile, checker, htedb);
		}
	}

	protected void executeCommand(String [] command, ExecutionContext exec, File lockFile) throws Exception {
		this.executeCommand(command, exec, null, lockFile, null, null, null, null, null);
	}
	
	protected void executeCommand(String [] command, ExecutionContext exec, File lockFile, String stdOutFile) throws Exception {
		this.executeCommand(command, exec, null, lockFile, stdOutFile, null, null, null, null);
	}
	
	protected void executeCommand(String [] command, ExecutionContext exec, File lockFile, String stdOutFile, String stdErrFile) throws Exception {
		this.executeCommand(command, exec, null, lockFile, stdOutFile, stdErrFile, null, null, null);
	}
	
	protected void executeCommand(String[] command, ExecutionContext exec,
			String[] environment, File lockFile, String stdOutFile,
			String stdErrFile, StringBuffer stdOut, StringBuffer stdErr,
			String StdInFile) throws Exception {
		
		use_hte = IBISKNIMENodesPlugin.getDefault().getHTEPreference();
		threshold_value = threshold.getIntValue();
		db_file = IBISKNIMENodesPlugin.getDefault().getDBFilePreference();
		boolean notify = IBISKNIMENodesPlugin.getDefault().getNotifyPreference();
		if(notify) {
			email = IBISKNIMENodesPlugin.getDefault().getEmailPreference();
		} else {
			email = null;
		}
		

		this.count = 0;
		//prepare klock
		host_name = InetAddress.getLocalHost().getHostName();
		lockCommand = "";
		for (String s : command) {
			lockCommand += s;
		}

		LOGGER.info(node_name+" is executed with: use_hte="+use_hte+" threshold="+threshold_value);
		
		if(lockFile == null) {
			lockFile = new File(defaultLockFile);
		}
		LOGGER.info("klock file can be found in "+lockFile);
		boolean terminationState = SuccessfulRunChecker
				.hasTerminatedSuccessfully(lockFile, lockCommand);
		LOGGER.info("Successful termination state: " + terminationState);

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
				htedb = new HTEDBHandler(db_file, LOGGER);
				if(!htedb.checkSchema()) throw new SQLException();
			} catch (SQLException e) {
				LOGGER.error("Connection to database could not be established: "+e.getMessage());
				use_hte = false;
			}
		}
		
		if (use_hte) {
			exec_id = htedb.insertNewHTExecution(lockCommand, node_name, host_name, threshold_value);
			if(stdErr == null) {
				stdErr = new StringBuffer();
			}
			recExecuteCommand(command, exec, environment, stdOutFile,
					stdErrFile, stdOut, stdErr, StdInFile, checker, htedb);
			htedb.closeConnection();
		} else {
			//HTE is not used
			LOGGER.info("HTE is not used");
			Executor.executeCommand(command, exec, environment, LOGGER,
					stdOutFile, stdErrFile, stdOut, stdErr, StdInFile);
			checker.writeOK();
			checker.finalize();
		}
	}
	
	private void sendMail(String content) {

		Properties properties = System.getProperties();
		properties.setProperty("mail.smtp.host", EMAILHOST);
		Session session = Session.getDefaultInstance(properties);

		try {
			MimeMessage message = new MimeMessage(session);
			message.setFrom(new InternetAddress(EMAILSENDER));
			message.addRecipient(Message.RecipientType.TO, new InternetAddress(
					email));
			message.setSubject(HEADER);
			message.setText(content);
			Transport.send(message);
		} catch (MessagingException mex) {
			mex.printStackTrace();
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

/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.InetAddress;
import java.sql.SQLException;
import java.util.Properties;

import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.ModelContent;
import org.knime.core.node.ModelContentRO;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.port.PortType;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.SettingsStorageNodeModel;
import de.helmholtz_muenchen.ibis.utils.threads.ExecuteThread;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;

/**
 * This is the model implementation of HTExecutorNode.
 * 
 * 
 * @author Tim Jeske
 */
public abstract class HTExecutorNodeModel extends SettingsStorageNodeModel {
	
	private static final NodeLogger LOGGER  = NodeLogger.getLogger(HTExecutorNodeModel.class);
	private StringBuffer HTEOUT 			= new StringBuffer("");
	private StringBuffer HTEERR 			= new StringBuffer("");
	
	//Internals
	private static final String INTERNAL_MODEL = "internalModel";
	private static final String INTERNALS_FILE = "hte_internals.xml";
	
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
	private String emailsender = "";
	private String emailhost = "";
	
	private final String HEADER = "IBIS KNIME Nodes Notification";
	
	//internal flags
	private boolean use_hte = false;
	
	static final String CFGKEY_DEFAULT_THRESHOLD = "threshold";
	static final int DEFAULT_THRESHOLD = 1;
	
	static final String CFGKEY_USE_PREF = "use_pref";
	private final SettingsModelBoolean m_use_pref = new SettingsModelBoolean(CFGKEY_USE_PREF, true);
	
	private final SettingsModelIntegerBounded threshold = new SettingsModelIntegerBounded(
			HTExecutorNodeModel.CFGKEY_DEFAULT_THRESHOLD, DEFAULT_THRESHOLD,1,Integer.MAX_VALUE);

	protected HTExecutorNodeModel(PortType[] inPortTypes,
			PortType[] outPortTypes) {
		super(inPortTypes, outPortTypes);
		init();
	}

	protected HTExecutorNodeModel(int nrInDataPorts, int nrOutDataPorts) {
		super(nrInDataPorts, nrOutDataPorts);
		init();
	}
	
	public void init() {
		addSetting(m_use_pref);
		addSetting(threshold);
	}
	
	@Override
	protected void reset(){
		resetView();
	}
	
	/**
	 * Resets Stdout/Stderr node views
	 */
	protected void resetView(){
		HTEOUT.setLength(0);
		HTEERR.setLength(0);
	}
	
	private void recExecuteCommand(String[] command, ExecutionContext exec,
			String[] environment, String stdOutFile,
			String stdErrFile, StringBuffer stdOut, StringBuffer stdErr,
			String StdInFile, SuccessfulRunChecker checker, HTEDBHandler htedb) throws Exception {

		if(stdErr==null){
			this.HTEERR = new StringBuffer();
			stdErr		= this.HTEERR;
		}else{
			this.HTEERR=stdErr;
		}
		
		int exitcode = 0;
		String err_msg = "";
		
		if (count < threshold_value) {
			count++;
			htedb.updateCount(exec_id, count);
			exitcode = Executor.executeCommandWithExitCode(command, exec,
					environment,LOGGER, stdOutFile, stdErrFile, stdOut,
					stdErr, StdInFile,HTEOUT);
			err_msg = stdErr.toString();
			
//			for(String c: command) {
//				if(c.contains(IBISKNIMENodesPlugin.GATK)) {
//					err_msg = parseGATKError(err_msg);
//					break;
//				}
//			}
			
			err_msg = err_msg.trim();

			stdErr = new StringBuffer();
			if (exitcode == 0) {
				htedb.writeSuccess(exec_id);
				checker.writeOK();
				checker.finalize();
				return;
			} else {
				writeHTELog(err_msg,"Error");
				
				htedb.writeError(exec_id, err_msg);
			}
		}

		if (count >= threshold_value) {
			String cmd = StringUtils.join(command, " ");
			
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
	
//	private String parseGATKError(String err) {
//		String result = "";
//		
//		for(String line: err.split(System.getProperty("line.separator"))) {
//			if(line.contains("ERROR MESSAGE")) {
//				result += line + System.getProperty("line.separator");
//			}
//		}
//		return result;
//	}

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
		
		if(stdErr==null){
			this.HTEERR = new StringBuffer();
			stdErr		= this.HTEERR;
		}else{
			this.HTEERR=stdErr;
		}
		
		use_hte = IBISKNIMENodesPlugin.getBooleanPreference(IBISKNIMENodesPlugin.USE_HTE);
		threshold_value = threshold.getIntValue();
		db_file = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.DB_FILE);
		boolean notify = IBISKNIMENodesPlugin.getBooleanPreference(IBISKNIMENodesPlugin.NOTIFY);
		if(notify) {
			email = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.EMAIL_RECEIVER);
			emailhost = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.EMAIL_HOST);
			emailsender = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.EMAIL_SENDER);
		} else {
			email = null;
		}
		

		this.count = 0;
		//prepare klock
		host_name = InetAddress.getLocalHost().getHostName();
		lockCommand = "";
		for (String s : command) {
			lockCommand = lockCommand + " " + s;
		}
		lockCommand = lockCommand.trim();

		writeHTELog(node_name+" is executed with: use_hte="+use_hte+" threshold="+threshold_value+" Host"+host_name,"Info");
		
		if(lockFile == null) {
			lockFile = new File(defaultLockFile);
		}
		
		writeHTELog("klock file can be found in "+lockFile,"Info");
		boolean terminationState = SuccessfulRunChecker
				.hasTerminatedSuccessfully(lockFile, lockCommand);
		writeHTELog("Successful termination state: " + terminationState,"Info");

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
				writeHTELog("Connection to database could not be established: "+e.getMessage(),"Error");
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
			writeHTELog("HTE is not used","Info");
			Executor.executeCommand(command, exec, environment, LOGGER,
					stdOutFile, stdErrFile, stdOut, stdErr, StdInFile,HTEOUT);
			checker.writeOK();
			checker.finalize();
		}
	}
	
	private void sendMail(String content) {

		Properties properties = System.getProperties();
		properties.setProperty("mail.smtp.host", emailhost);
		Session session = Session.getDefaultInstance(properties);

		try {
			MimeMessage message = new MimeMessage(session);
			message.setFrom(new InternetAddress(emailsender));
			message.addRecipient(Message.RecipientType.TO, new InternetAddress(
					email));
			message.setSubject(HEADER);
			message.setText(content);
			Transport.send(message);
		} catch (MessagingException mex) {
			mex.printStackTrace();
		}
	}

	public String getHTEOUT() {
		return HTEOUT.toString();
	}
	public String getHTEERR() {
		return HTEERR.toString();
	}
	
	private void writeHTELog(String Log, String Type){
		if(Type.equals("Info")){
			LOGGER.info(Log);
		}else if(Type.equals("Error")){
			LOGGER.error(Log);
		}else{
			throw new InternalError("Incorrect call of writeHTELog. Unknown Log Type: "+Type);
		}
		
		HTEOUT.append(Log+"\n");
	}
	
	
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    	 File file = new File(internDir, INTERNALS_FILE);
         FileInputStream fis = new FileInputStream(file);
         ModelContentRO modelContent = ModelContent.loadFromXML(fis);
         try {
			HTEOUT.append(modelContent.getString("ViewSTDOUT"));
			HTEERR.append(modelContent.getString("ViewSTDERR"));
		} catch (InvalidSettingsException e) {
			e.printStackTrace();
		}
    }
    

    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    	ModelContent modelContent = new ModelContent(INTERNAL_MODEL);
    	modelContent.addString("ViewSTDOUT", getHTEOUT());
    	modelContent.addString("ViewSTDERR", getHTEERR());
    	File file = new File(internDir, INTERNALS_FILE);
    	FileOutputStream fos = new FileOutputStream(file);
        modelContent.saveToXML(fos);
    }
	
	
	
//	/**
//	 * {@inheritDoc}
//	 */
//	@Override
//	protected void saveSettingsTo(final NodeSettingsWO settings) {
//		threshold.saveSettingsTo(settings);
//		m_use_pref.saveSettingsTo(settings);
//	}
//
//	/**
//	 * {@inheritDoc}
//	 */
//	@Override
//	protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
//			throws InvalidSettingsException {
//		threshold.loadSettingsFrom(settings);
//		m_use_pref.loadSettingsFrom(settings);
//	}
//
//	/**
//	 * {@inheritDoc}
//	 */
//	@Override
//	protected void validateSettings(final NodeSettingsRO settings)
//			throws InvalidSettingsException {
//		threshold.validateSettings(settings);
//		m_use_pref.validateSettings(settings);
//	}
}

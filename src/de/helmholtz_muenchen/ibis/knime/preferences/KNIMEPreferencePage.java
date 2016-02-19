
package de.helmholtz_muenchen.ibis.knime.preferences;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.sql.SQLException;
import java.util.HashMap;

import javax.mail.internet.AddressException;
import javax.mail.internet.InternetAddress;
import javax.swing.JOptionPane;

import org.apache.commons.io.FileUtils;
import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.jface.preference.PreferencePage;
import org.eclipse.swt.SWT;
import org.eclipse.swt.events.SelectionAdapter;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.layout.GridData;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Control;
import org.eclipse.swt.widgets.DirectoryDialog;
import org.eclipse.swt.widgets.Event;
import org.eclipse.swt.widgets.FileDialog;
import org.eclipse.swt.widgets.Group;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Listener;
import org.eclipse.swt.widgets.Shell;
import org.eclipse.swt.widgets.Table;
import org.eclipse.swt.widgets.TableColumn;
import org.eclipse.swt.widgets.TableItem;
import org.eclipse.swt.widgets.Text;
import org.eclipse.ui.IWorkbench;
import org.eclipse.ui.IWorkbenchPreferencePage;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.BinaryHandler;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTEDBHandler;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;

/**
 * @author hastreiter
 * @author jeske
 */
public class KNIMEPreferencePage extends PreferencePage implements
        IWorkbenchPreferencePage {

	public static final String GATK_ID = "GenomeAnalysisTK.jar";
	private static final String DOWNLOAD_PATH = "ftp://ftpmips.helmholtz-muenchen.de/Incoming/KNIME_BIN/";
	public static final HashMap<String, Boolean> TOOLS;
	static {
		TOOLS = new HashMap<>();
		TOOLS.put("bfast",true);
		TOOLS.put("bowtie2",true);
		TOOLS.put("bwa",true);
		TOOLS.put("featureCounts",true);
		TOOLS.put("pindel",true);
		TOOLS.put("pindel2vcf",true);
		TOOLS.put("samtools",true);
		TOOLS.put("segemehl.x",true);
		TOOLS.put("STAR",true);
		TOOLS.put(GATK_ID,false);
		TOOLS.put("variant_effect_predictor.pl",false);
		TOOLS.put("filter_vep.pl", false);
		TOOLS.put("vcftools",false);
	}
	
	public static String REF_GENOME;
	
	public static boolean USE_HTE;
	public static String THRESHOLD;
	public static String DB_FILE;
	public static boolean NOTIFY;
	public static String EMAILSENDER;
	public static String EMAILHOST;
	public static String EMAILRECEIVER;
	
	private Text refGenome;
	private Text thresholdText;
	private Text dbFile;
	private Text email_sender, email_host, email_receiver;
	
	private Button checkHTE;
	private Button checkNotify;
	
	private Table table;

	private static final NodeLogger LOGGER = NodeLogger.getLogger(KNIMEPreferencePage.class);

	public KNIMEPreferencePage() {
		super();

		// Set the preference store for the preference page.
		IPreferenceStore store =
				IBISKNIMENodesPlugin.getDefault().getPreferenceStore();
		setPreferenceStore(store);
	}

	/**
	 * @see org.eclipse.jface.preference.
	 * PreferencePage#createContents(Composite)
	 */
	protected Control createContents(Composite parent) {
		
		Composite top = new Composite(parent, SWT.LEFT);
		top.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		top.setLayout(new GridLayout());
		
		//bin preferences
		Group binPrefs = new Group(top, SWT.NONE);
		binPrefs.setText("Binary preferences");
		binPrefs.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		binPrefs.setLayout(new GridLayout());

		GridData gd_table = new GridData(SWT.FILL, SWT.TOP, true, false);
		gd_table.heightHint = 200;
		table = new Table(binPrefs, SWT.BORDER | SWT.V_SCROLL | SWT.H_SCROLL);
		table.setLayoutData(gd_table);
		table.setHeaderVisible(true);
		String [] titles = {"Tool", "Path to binary"};
		
		for(int i = 0; i < titles.length; i++) {
			TableColumn col = new TableColumn(table, SWT.NULL);
		
			col.setText(titles[i]);
		}
		
		for(String key: TOOLS.keySet()) {
			TableItem item = new TableItem(table, SWT.NULL);
			item.setText(0,key);
			item.setText(1,IBISKNIMENodesPlugin.getDefault().getToolPathPreference(key));
		}
		
		for(int i = 0; i < titles.length; i++) {
			table.getColumn(i).pack();
		}

		
		table.addListener(SWT.Selection, new Listener() {
			public void handleEvent(Event event) {
				LOGGER.debug("You selected " + event.item);
				}
		});
		
		Composite searchDownloadEdit = new Composite(binPrefs, SWT.NONE);
		GridLayout downloadLayout = new GridLayout();
		downloadLayout.numColumns = 3;
		searchDownloadEdit.setLayout(downloadLayout);
		
		Button browseSearchDir = new Button(searchDownloadEdit, SWT.NONE);
		browseSearchDir.setText("Search in directory");
		final Shell shell = new Shell(parent.getDisplay());
		browseSearchDir.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				selectSearchDir(shell);
			}
		});
		
		Button downloader = new Button(searchDownloadEdit, SWT.NONE);
		downloader.setText("Download missing binaries");
		final Shell shell1 = new Shell(parent.getDisplay());
		downloader.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				downloadBinaries(shell1);
			}
		});
		
		Button edit = new Button(searchDownloadEdit, SWT.NONE);
		edit.setText("Edit binary");
		final Shell shell2 = new Shell(parent.getDisplay());
		edit.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				editBinary(shell2);
			}
		});
		
		//Reference genome
		Group ref_genome = new Group(top,SWT.NONE);
		ref_genome.setText("Reference (genome) sequence");
		ref_genome.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		GridLayout refLayout = new GridLayout();
		refLayout.numColumns = 3;
		ref_genome.setLayout(refLayout);
		
		REF_GENOME = IBISKNIMENodesPlugin.getDefault().getRefGenomePreference();
		
		Label refGenomeLabel = new Label(ref_genome,SWT.LEFT);
		refGenomeLabel.setText("File path:");
		
		refGenome = new Text(ref_genome,SWT.BORDER);
		refGenome.setText(REF_GENOME);
		refGenome.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		
		Button browseRefGenome = new Button(ref_genome, SWT.RIGHT);
		browseRefGenome.setText("Browse");
		final Shell shell5 = new Shell(parent.getDisplay());
		browseRefGenome.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				selectRefGenome(shell5);
			}
		});
		
		//HTE preferences
		Group htePrefs = new Group(top,SWT.NONE);
		htePrefs.setText("High throughput execution preferences");
		htePrefs.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		htePrefs.setLayout(new GridLayout());
		
		Composite use_hte = new Composite(htePrefs,SWT.LEFT);
		use_hte.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		GridLayout hteLayout = new GridLayout();
		hteLayout.numColumns = 3;
		use_hte.setLayout(hteLayout);
		
		USE_HTE = IBISKNIMENodesPlugin.getDefault().getHTEPreference();
		checkHTE = new Button(use_hte,SWT.CHECK);
		checkHTE.setText("Use HTE?");
		checkHTE.setSelection(USE_HTE);
		
		
		Label thresholdLabel = new Label(use_hte,SWT.RIGHT);
		thresholdLabel.setText("Global threshold:");
		
		THRESHOLD = IBISKNIMENodesPlugin.getDefault().getThresholdPreference();

		thresholdText = new Text(use_hte, SWT.BORDER);
		thresholdText.setText(THRESHOLD);
		thresholdText.setEnabled(USE_HTE);
		
		Label dbFileLabel = new Label(use_hte, SWT.LEFT);
		dbFileLabel.setText("Use existing db file:");
		
		DB_FILE = IBISKNIMENodesPlugin.getDefault().getDBFilePreference();
		
		dbFile = new Text(use_hte,SWT.BORDER);
		dbFile.setText(DB_FILE);
		dbFile.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		dbFile.setEnabled(USE_HTE);
		dbFile.setEditable(false);
		
		Button browseDBFile = new Button(use_hte, SWT.RIGHT);
		browseDBFile.setText("Browse");
		browseDBFile.setEnabled(USE_HTE);
		final Shell shell3 = new Shell(parent.getDisplay());
		browseDBFile.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				selectDBFile(shell3);
			}
		});
		
		Button createDBFile = new Button(use_hte, SWT.NONE);
		createDBFile.setText("Create new db file");
		createDBFile.setEnabled(USE_HTE);
		final Shell shell4 = new Shell(parent.getDisplay());
		createDBFile.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				createDBFile(shell4);
			}
		});
		
		checkHTE.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				Button btn = (Button) e.getSource();
				USE_HTE = btn.getSelection();
				thresholdText.setEnabled(USE_HTE);
				dbFile.setEnabled(USE_HTE);
				checkNotify.setEnabled(USE_HTE);
				checkNotify.setSelection(false);
				browseDBFile.setEnabled(USE_HTE);
				createDBFile.setEnabled(USE_HTE);
				NOTIFY = false;
				email_receiver.setEnabled(NOTIFY);
				email_host.setEnabled(NOTIFY);
				email_sender.setEnabled(NOTIFY);
			}
		});
		
		//email preferences
		Composite email = new Composite(htePrefs,SWT.LEFT);
		email.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		GridLayout emailLayout = new GridLayout();
		emailLayout.numColumns = 1;
		email.setLayout(emailLayout);
		
		NOTIFY = IBISKNIMENodesPlugin.getDefault().getNotifyPreference();
		checkNotify = new Button(email,SWT.CHECK);
		checkNotify.setText("Email notification?");
		checkNotify.setSelection(NOTIFY);
		checkNotify.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				Button btn = (Button) e.getSource();
				NOTIFY = btn.getSelection();
				email_receiver.setEnabled(NOTIFY);
				email_host.setEnabled(NOTIFY);
				email_sender.setEnabled(NOTIFY);
			}
		});
		
		EMAILHOST = IBISKNIMENodesPlugin.getDefault().getEmailPreference(IBISKNIMENodesPlugin.EMAIL_HOST);
		
		Label emailHost = new Label(email, SWT.LEFT);
		emailHost.setText("Email host:");
		email_host = new Text(email,SWT.BORDER);
		email_host.setText(EMAILHOST);
		email_host.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		email_host.setEnabled(NOTIFY);
		
		EMAILSENDER = IBISKNIMENodesPlugin.getDefault().getEmailPreference(IBISKNIMENodesPlugin.EMAIL_SENDER);
		
		Label emailSender = new Label(email, SWT.LEFT);
		emailSender.setText("Email sender:");
		email_sender = new Text(email,SWT.BORDER);
		email_sender.setText(EMAILSENDER);
		email_sender.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		email_sender.setEnabled(NOTIFY);
		
		EMAILRECEIVER = IBISKNIMENodesPlugin.getDefault().getEmailPreference(IBISKNIMENodesPlugin.EMAIL_RECEIVER);
		
		Label emailReceiver = new Label(email, SWT.LEFT);
		emailReceiver.setText("Email address (send to):");
		email_receiver = new Text(email,SWT.BORDER);
		email_receiver.setText(EMAILRECEIVER);
		email_receiver.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		email_receiver.setEnabled(NOTIFY);

		return top;
	}

	/**
	 * @see IWorkbenchPreferencePage#init
	 */	
	public void init(IWorkbench wb) {	
	}
	
	/*
	 * The user has pressed "Restore defaults".
	 * Restore all default preferences.
	 */
	protected void performDefaults() {
		super.performDefaults();
		
		IBISKNIMENodesPlugin iknp = IBISKNIMENodesPlugin.getDefault();
		
		for(String s: TOOLS.keySet()) {
			iknp.setToolPathPreference(s, "");
		}
		
		for(TableItem i: table.getItems()) {
			i.setText(1,"");
		}
		
		REF_GENOME = IBISKNIMENodesPlugin.REF_GENOME_DEFAULT;
		refGenome.setText(REF_GENOME);
		iknp.setRefGenomePreference(REF_GENOME);
		
		USE_HTE = IBISKNIMENodesPlugin.HTE_DEFAULT;
		checkHTE.setSelection(USE_HTE);
		iknp.setHTEPreference(USE_HTE);
		
		THRESHOLD = IBISKNIMENodesPlugin.THRESHOLD_DEFAULT+"";
		thresholdText.setText(THRESHOLD);
		iknp.setThresholdPreference(THRESHOLD);
		
		DB_FILE = IBISKNIMENodesPlugin.DB_FILE_DEFAULT;
		dbFile.setText(DB_FILE);
		iknp.setDBFilePreference(DB_FILE);
		
		NOTIFY = IBISKNIMENodesPlugin.NOTIFY_DEFAULT;
		checkNotify.setSelection(NOTIFY);
		iknp.setNotifyPreference(NOTIFY);
		
		EMAILHOST = IBISKNIMENodesPlugin.EMAIL_HOST_DEFAULT;
		email_host.setText(EMAILHOST);
		iknp.setEmailPreference(IBISKNIMENodesPlugin.EMAIL_HOST, EMAILHOST);
		
		EMAILSENDER = IBISKNIMENodesPlugin.EMAIL_SENDER_DEFAULT;
		email_sender.setText(EMAILSENDER);
		iknp.setEmailPreference(IBISKNIMENodesPlugin.EMAIL_SENDER, EMAILSENDER);
		
		EMAILRECEIVER = IBISKNIMENodesPlugin.EMAIL_RECEIVER_DEFAULT;
		email_receiver.setText(EMAILRECEIVER);
		iknp.setEmailPreference(IBISKNIMENodesPlugin.EMAIL_RECEIVER, EMAILRECEIVER);
	}
	
	/*
	 * The user has pressed Ok or Apply. Store/apply 
	 * this page's values appropriately.
	 */	
	public boolean performOk() {
		
		IBISKNIMENodesPlugin iknp = IBISKNIMENodesPlugin.getDefault();
		
		REF_GENOME = refGenome.getText();
		if(!REF_GENOME.equals("") && !FileValidator.checkFastaFormat(REF_GENOME)) {
			JOptionPane.showMessageDialog(null,
					"Reference (genome) sequence file is not in FastA format or does not contain nucleotide sequences!",
				    "Error",
				    JOptionPane.ERROR_MESSAGE);
			return false;
		}
		iknp.setRefGenomePreference(REF_GENOME);
		LOGGER.debug("Setting REF_GENOME to: "+REF_GENOME);
		
		iknp.setHTEPreference(USE_HTE);
		LOGGER.debug("Setting USE_HTE to: "+USE_HTE);
		
		if(!USE_HTE) {
			iknp.setNotifyPreference(false);
			return super.performOk();
		}
		
		THRESHOLD = thresholdText.getText();
		try{
			int n = Integer.parseInt(THRESHOLD);
			if(n<0) {
				throw(new NumberFormatException());
			}
		} catch (NumberFormatException e){
			JOptionPane.showMessageDialog(null,
				    "Threshold has to be a positive integer.",
				    "Error",
				    JOptionPane.ERROR_MESSAGE);
			return false;
		}
		iknp.setThresholdPreference(THRESHOLD);
		LOGGER.debug("Setting THREHOLD to: "+THRESHOLD);
		
		if(DB_FILE.equals("")) {
			JOptionPane.showMessageDialog(null,
				    "HTE requires a SQLite database. Please create a new one or choose an existing one.",
				    "Error",
				    JOptionPane.ERROR_MESSAGE);
			return false;
		}
		iknp.setDBFilePreference(DB_FILE);
		LOGGER.debug("Setting DB_FILE to: "+DB_FILE);
		
		iknp.setNotifyPreference(NOTIFY);
		
		EMAILHOST = email_host.getText();
		iknp.setEmailPreference(IBISKNIMENodesPlugin.EMAIL_HOST, EMAILHOST);
		
		EMAILSENDER = email_sender.getText();
		iknp.setEmailPreference(IBISKNIMENodesPlugin.EMAIL_SENDER, EMAILSENDER);
		
		EMAILRECEIVER = email_receiver.getText();
		iknp.setEmailPreference(IBISKNIMENodesPlugin.EMAIL_RECEIVER, EMAILRECEIVER);
		
		if(NOTIFY) {
			
			if(EMAILHOST.equals("")) {
				JOptionPane.showMessageDialog(null,
					    "Email host is required for Email notification!",
					    "Error",
					    JOptionPane.ERROR_MESSAGE);
				return false;
			}
			
			if(EMAILSENDER.equals("")) {
				JOptionPane.showMessageDialog(null,
					    "Email sender is required for Email notification!",
					    "Error",
					    JOptionPane.ERROR_MESSAGE);
				return false;
			}
			
			try {
				new InternetAddress(EMAILRECEIVER).validate();
			} catch (AddressException e) {
				JOptionPane.showMessageDialog(null,
						"Enter valid Email address.",
						"Error",
						JOptionPane.ERROR_MESSAGE);
				return false;
			}
		}
		
		return super.performOk();
	}
	
	private void downloadBinaries(Shell shell){
		
		IBISKNIMENodesPlugin iknp = IBISKNIMENodesPlugin.getDefault();
		
		DirectoryDialog dlg = new DirectoryDialog(shell);
		dlg.setText("Choose directory in which tool binaries will be stored");
		dlg.setFilterPath("~/");
		String dir = dlg.open();
		
		//do nothing if DirectoryDialog was cancelled
		if(dir==null) return;
		
		try {
			CheckUtils.checkDestinationDirectory(dir);
		} catch (InvalidSettingsException e1) {
			JOptionPane.showMessageDialog(null,
				    "This directory cannot be used: "+e1.getMessage(),
				    "Error",
				    JOptionPane.ERROR_MESSAGE);
			return;
		}
		
			
		for (String tool : TOOLS.keySet()) {
			if (TOOLS.get(tool)) {
				if (iknp.getToolPathPreference(tool).equals("")) {
					try {
						File f = new File(dir + "/" + tool);
						if (!f.exists()) {
							FileUtils.copyURLToFile(new URL(DOWNLOAD_PATH + tool), f);
							f.setExecutable(true, false);
							iknp.setToolPathPreference(tool, dir + "/"+ tool);
						}
					} catch (IOException e) {
						LOGGER.error("Downloading " + tool+ " failed! Message: " + e.getMessage());
					}
				}
			}
		}
		
		for(TableItem i: table.getItems()) {
			i.setText(1,iknp.getToolPathPreference(i.getText(0)));
		}
	}
	
	private void selectSearchDir(Shell shell){
		
		IBISKNIMENodesPlugin iknp = IBISKNIMENodesPlugin.getDefault();
		
		DirectoryDialog dirlg = new DirectoryDialog(shell);
		dirlg.setText("Choose directory in which tool binaries shall be searched");
		dirlg.setFilterPath("~/");
		String dir = dirlg.open();
		
		//do nothing if DirectoryDialog was cancelled
		if(dir==null) return;

	    Thread t = new Thread(new Runnable() {
	      public void run() {
	    	  JOptionPane.showMessageDialog(null,
	  			    "Please wait while directories are searched. The file dialog will be closed automatically.",
	  			    "Information",
	  			    JOptionPane.INFORMATION_MESSAGE);
	      }
	    });
	    t.start();
		
		for(String s: TOOLS.keySet()) {
			if(iknp.getToolPathPreference(s).equals("")) {
				String path = BinaryHandler.checkToolAvailability(s, dir);
				if(path != null) {
					try {
						CheckUtils.checkSourceFile(path);
						iknp.setToolPathPreference(s, path);
					} catch (InvalidSettingsException e) {
						iknp.setToolPathPreference(s, "");
					}
					
				}
			}
		}
		
		
		for(TableItem i: table.getItems()) {
			i.setText(1,iknp.getToolPathPreference(i.getText(0)));
		}
	}
	
	private void editBinary(Shell shell) {
		
		IBISKNIMENodesPlugin iknp = IBISKNIMENodesPlugin.getDefault();
		
		if(table.getSelectionCount()==0) {
			return;
		}
		TableItem item = table.getSelection()[0];
		
		FileDialog dlg = new FileDialog(shell);
		dlg.setText("Select directory path to "+item.getText(0));
		String path = System.getProperty("user.home");
		dlg.setFilterPath(path);
		String file = dlg.open();
		
		//do nothing if FileDialog was cancelled
		if(file==null) return;
		
		try {
			CheckUtils.checkSourceFile(file);
			iknp.setToolPathPreference(item.getText(0), file);
			item.setText(1,file);
		} catch (InvalidSettingsException e) {
			JOptionPane.showMessageDialog(null,
				    "This file cannot be used: "+e.getMessage(),
				    "Error",
				    JOptionPane.ERROR_MESSAGE);
			return;
		}
	}
	
	private void selectRefGenome(Shell shell) {
		FileDialog fdl = new FileDialog(shell);
		fdl.setText("Select reference genome");
		fdl.setFilterExtensions(new String[]{"*.fa;*.fasta"});
		fdl.setFilterPath("~/");
		String path = fdl.open();
		
		//do nothing if FileDialog was cancelled
		if(path==null) return;
		
		REF_GENOME = path;
		refGenome.setText(REF_GENOME);
	}
	
	private void selectDBFile(Shell shell) {
		FileDialog fdl = new FileDialog(shell);
		fdl.setText("Selecte HTE database file");
		fdl.setFilterExtensions(new String[]{"*.db","*.sql"});
		fdl.setFilterPath("~/");
		String path = fdl.open();
		
		//do nothing if FileDialog was cancelled
		if(path==null) return;
		
		HTEDBHandler htedb;
		try {
			htedb = new HTEDBHandler(path, null);
			if(htedb.checkSchema()) {
				DB_FILE = path;
			} else {
				JOptionPane.showMessageDialog(null,
					    "Database schema of selected database is not valid.",
					    "Error",
					    JOptionPane.ERROR_MESSAGE);
				DB_FILE = "";
			}
			htedb.closeConnection();
			dbFile.setText(DB_FILE);
		} catch (SQLException e) {
			JOptionPane.showMessageDialog(null,
					"Connecting to DB failed!"+ System.getProperty("line.separator")+ "Error message: "+e.getMessage(),
					"Error",
					JOptionPane.ERROR_MESSAGE);
			DB_FILE="";
			dbFile.setText(DB_FILE);
		}
	}
	
	public void createDBFile(Shell shell) {
		DirectoryDialog dlg = new DirectoryDialog(shell);
		dlg.setText("Choose directory in which database file will be stored");
		dlg.setFilterPath("~/");
		String dir = dlg.open();
		
		//do nothing if DirectoryDialog was cancelled
		if(dir==null) return;
		
		try {
			HTEDBHandler htedb = new HTEDBHandler(dir+"/hte.db",null);
			htedb.createDB();
			htedb.closeConnection();
			DB_FILE = dir+"/hte.db";
		} catch (SQLException e) {
			JOptionPane.showMessageDialog(null,
					"Creating DB failed!"+ System.getProperty("line.separator")+ "Error message: "+e.getMessage(),
					"Error",
					JOptionPane.ERROR_MESSAGE);
			DB_FILE="";
			
		}
		dbFile.setText(DB_FILE);	
	}
}

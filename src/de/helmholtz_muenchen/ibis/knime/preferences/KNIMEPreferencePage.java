
package de.helmholtz_muenchen.ibis.knime.preferences;


import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.sql.SQLException;

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
import org.eclipse.swt.widgets.FileDialog;
import org.eclipse.swt.widgets.Group;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Shell;
import org.eclipse.swt.widgets.Text;
import org.eclipse.ui.IWorkbench;
import org.eclipse.ui.IWorkbenchPreferencePage;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTEDBHandler;

/**
 * @author hastreiter
 */
public class KNIMEPreferencePage extends PreferencePage implements
        IWorkbenchPreferencePage {

	private final String DOWNLOAD_PATH = "ftp://ftpmips.helmholtz-muenchen.de/Incoming/KNIME_BIN/";
	private final String [] TOOLS = {"bwa","pindel","pindel2vcf","GenomeAnalysisTK.jar"};
	
	public static String TOOL_LOCATION;
	public static boolean USE_HTE;
	public static String THRESHOLD;
	public static String DB_FILE;
	public static boolean NOTIFY;
	public static String EMAIL;
	
	private Text binsDirectory;
	private Text thresholdText;
	private Text dbFile;
	private Text email;
	
	private Button checkHTE;
	private Button checkNotify;


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
		Group downloadBins = new Group(top, SWT.NONE);
		downloadBins.setText("Binary preferences");
		downloadBins.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		GridLayout downloadLayout = new GridLayout();
		downloadLayout.numColumns = 3;
		downloadBins.setLayout(downloadLayout);

		Label binsLabel = new Label(downloadBins, SWT.NONE);
		binsLabel.setText("Directory for tool binaries:");
		
		TOOL_LOCATION = IBISKNIMENodesPlugin.getDefault().getToolDirPreference();
		
		binsDirectory = new Text(downloadBins, SWT.BORDER);
		binsDirectory.setText(TOOL_LOCATION);
		binsDirectory.setEditable(false);
		
		Button browseBinDir = new Button(downloadBins, SWT.NONE);
		browseBinDir.setText("Browse");
		final Shell shell = new Shell(parent.getDisplay());
		browseBinDir.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				selectToolFolder(shell);
			}
		});

		Button downloader = new Button(downloadBins, SWT.NONE);
		downloader.setText("Download missing binaries");
		downloader.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				downloadBinaries();
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
		final Shell shell2 = new Shell(parent.getDisplay());
		browseDBFile.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				selectDBFile(shell2);
			}
		});
		
		Button createDBFile = new Button(use_hte, SWT.NONE);
		createDBFile.setText("Create new db file");
		final Shell shell3 = new Shell(parent.getDisplay());
		createDBFile.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				createDBFile(shell3);
			}
		});
		
		checkHTE.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				Button btn = (Button) e.getSource();
				USE_HTE = btn.getSelection();
				thresholdText.setEnabled(USE_HTE);
				dbFile.setEnabled(USE_HTE);
			}
		});
		
		Composite address = new Composite(htePrefs,SWT.LEFT);
		address.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		GridLayout emailLayout = new GridLayout();
		emailLayout.numColumns = 1;
		address.setLayout(emailLayout);
		
		NOTIFY = IBISKNIMENodesPlugin.getDefault().getNotifyPreference();
		checkNotify = new Button(address,SWT.CHECK);
		checkNotify.setText("Email notification?");
		checkNotify.setSelection(NOTIFY);
		checkNotify.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				Button btn = (Button) e.getSource();
				NOTIFY = btn.getSelection();
				email.setEnabled(NOTIFY);
			}
		});
		
		Label emailLabel = new Label(address, SWT.LEFT);
		emailLabel.setText("Email address:");
		
		EMAIL = IBISKNIMENodesPlugin.getDefault().getEmailPreference();
		
		email = new Text(address,SWT.BORDER);
		email.setText(EMAIL);
		email.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		email.setEnabled(NOTIFY);

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
		
		TOOL_LOCATION = IBISKNIMENodesPlugin.TOOL_DIR_DEFAULT;
		binsDirectory.setText(TOOL_LOCATION);
		iknp.setToolDirPreference(TOOL_LOCATION);
		
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
		
		EMAIL = IBISKNIMENodesPlugin.EMAIL_DEFAULT;
		email.setText(EMAIL);
		iknp.setEmailPreference(EMAIL);
		
	}
	
	/*
	 * The user has pressed Ok or Apply. Store/apply 
	 * this page's values appropriately.
	 */	
	public boolean performOk() {
		
		IBISKNIMENodesPlugin iknp = IBISKNIMENodesPlugin.getDefault();
		
		iknp.setToolDirPreference(TOOL_LOCATION);
		System.out.println("Setting TOOL_LOCATION to: "+TOOL_LOCATION);
		
		iknp.setHTEPreference(USE_HTE);
		System.out.println("Setting USE_HTE to: "+USE_HTE);
		
		if(!USE_HTE) {
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
		System.out.println("Setting THREHOLD to: "+THRESHOLD);
		
		if(DB_FILE.equals("")) {
			JOptionPane.showMessageDialog(null,
				    "Select valid SQLite database or create new one.",
				    "Error",
				    JOptionPane.ERROR_MESSAGE);
			return false;
		}
		iknp.setDBFilePreference(DB_FILE);
		System.out.println("Setting DB_FILE to: "+DB_FILE);
		
		iknp.setNotifyPreference(NOTIFY);
		
		EMAIL = email.getText();
		iknp.setEmailPreference(EMAIL);
		System.out.println("Setting EMAIL to: "+EMAIL);
		
		if(NOTIFY) {
			try {
				new InternetAddress(EMAIL).validate();
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
	
	private void downloadBinaries(){
		
		String dir = TOOL_LOCATION;
		System.out.println(dir);

		if (dir != null) {
			for(String t:TOOLS) {
				try {
					File f = new File(dir+"/"+t);
					if(!f.exists()) {
						FileUtils.copyURLToFile(new URL(DOWNLOAD_PATH+t), f);
						f.setExecutable(true,false);
					}
				} catch (IOException e) {
					JOptionPane.showMessageDialog(null,
							"Downloading "+t +" failed!"+ System.getProperty("line.separator")+ "Error message: "+e.getMessage(),
							"Error",
							JOptionPane.ERROR_MESSAGE);
					TOOL_LOCATION="";
					binsDirectory.setText(TOOL_LOCATION);
					break;
				}
			}
		}
	}
	
	private void selectToolFolder(Shell shell){
		
		DirectoryDialog dlg = new DirectoryDialog(shell);
		dlg.setText("Choose directory in which tool binaries will be stored");
		dlg.setFilterPath("~/");
		String dir = dlg.open();
		TOOL_LOCATION=dir;
		binsDirectory.setText(TOOL_LOCATION);
	}
	
	private void selectDBFile(Shell shell) {
		FileDialog fdl = new FileDialog(shell);
		fdl.setText("Selecte HTE database file");
		fdl.setFilterExtensions(new String[]{"*.db","*.sql"});
		fdl.setFilterPath("~/");
		String path = fdl.open();
		
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

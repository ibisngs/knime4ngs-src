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
package de.helmholtz_muenchen.ibis.knime.preferences;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.HashSet;

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
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTEDBHandler;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;


/**
 * @author Maximilian Hastreiter
 * @author Tim Jeske
 */
public class KNIMEPreferencePage extends PreferencePage implements
        IWorkbenchPreferencePage {

	
	private static final String DOWNLOAD_PATH = "ftp://ftpmips.helmholtz-muenchen.de/knime/";

	
	//define all dependencies for each tool (recursive dependencies are not checked!)
	public static final HashMap<String, HashSet<String>> DEPENDENCIES;
	static {
		DEPENDENCIES = new HashMap<>();
		HashSet<String> bowtie2_dep = new HashSet<>();
		bowtie2_dep.add(IBISKNIMENodesPlugin.BOWTIE2_BUILD);
		bowtie2_dep.add(IBISKNIMENodesPlugin.BOWTIE2_BUILD_L);
		bowtie2_dep.add(IBISKNIMENodesPlugin.BOWTIE2_BUILD_S);
		bowtie2_dep.add(IBISKNIMENodesPlugin.BOWTIE2_ALIGN_L);
		bowtie2_dep.add(IBISKNIMENodesPlugin.BOWTIE2_ALIGN_S);
		DEPENDENCIES.put(IBISKNIMENodesPlugin.BOWTIE2, bowtie2_dep);
	}
	
	public static String REF_GENOME, RES_HAPMAP, RES_OMNI, RES_1000G_SNPS, RES_1000G_INDELS, RES_DBSNP, RES_MILLS;
	public static boolean OVERWRITE, USE_HTE;
	public static String THRESHOLD;
	public static String DB_FILE;
	public static boolean NOTIFY;
	public static String EMAILSENDER, EMAILHOST, EMAILRECEIVER;
	
	private Text refGenome, res_hapmap, res_omni, res_1000G_SNPS, res_1000G_Indels, res_dbsnp, res_mills; 
	private Text thresholdText;
	private Text dbFile;
	private Text email_sender, email_host, email_receiver;
	
	private Button checkOverwrite;
	private Button checkHTE;
	private Button checkNotify;
	
	public Table table;
	
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
		
		final Shell shell = new Shell(parent.getDisplay());
		
		Composite top = new Composite(parent, SWT.LEFT);
		top.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		top.setLayout(new GridLayout());
		
		//overwrite preferences
		Group outPrefs = new Group(top,SWT.NONE);
		outPrefs.setText("Output preferences");
		outPrefs.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		outPrefs.setLayout(new GridLayout());
		
		Composite overwrite = new Composite(outPrefs,SWT.LEFT);
		overwrite.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		GridLayout overwriteLayout = new GridLayout();
		overwriteLayout.numColumns = 1;
		overwrite.setLayout(overwriteLayout);
		
		OVERWRITE = IBISKNIMENodesPlugin.getBooleanPreference(IBISKNIMENodesPlugin.OVERWRITE);
		checkOverwrite = new Button(overwrite,SWT.CHECK);
		checkOverwrite.setText("Overwrite existing files?");
		checkOverwrite.setSelection(OVERWRITE);
		
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
		
		for(String key: IBISKNIMENodesPlugin.TOOLS.keySet()) {
			TableItem item = new TableItem(table, SWT.NULL);
			item.setText(0,key);
			item.setText(1,IBISKNIMENodesPlugin.getStringPreference(key));
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
		downloadLayout.numColumns = 4;
		searchDownloadEdit.setLayout(downloadLayout);
		
		Button browseSearchDir = new Button(searchDownloadEdit, SWT.NONE);
		browseSearchDir.setText("Search in directory");
		browseSearchDir.setEnabled(!IBISKNIMENodesPlugin.getDefault().isSearching());
		
		Button cancelSearch = new Button(searchDownloadEdit, SWT.NONE);
		cancelSearch.setText("Cancel search");
		cancelSearch.setEnabled(IBISKNIMENodesPlugin.getDefault().isSearching());
		
		Button downloader = new Button(searchDownloadEdit, SWT.NONE);
		downloader.setText("Download missing binaries");

		Button edit = new Button(searchDownloadEdit, SWT.NONE);
		edit.setText("Edit binary");
		
		//Reference genome and resource files
		GridLayout refLayout = new GridLayout();
		refLayout.numColumns = 3;
		
		Group ref_res_group = new Group(top,SWT.NONE);
		ref_res_group.setText("Reference sequence and resource data sets");
		ref_res_group.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		ref_res_group.setLayout(refLayout);
		
		REF_GENOME = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.REF_GENOME);
		
		Label refGenomeLabel = new Label(ref_res_group,SWT.LEFT);
		refGenomeLabel.setText("Reference sequence:");
		
		refGenome = new Text(ref_res_group,SWT.BORDER);
		refGenome.setText(REF_GENOME);
		refGenome.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		
		Button browseRefGenome = new Button(ref_res_group, SWT.RIGHT);
		browseRefGenome.setText("Browse");
		
		RES_HAPMAP = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.RES_HAPMAP);
		Label resHapmapLabel = new Label(ref_res_group,SWT.LEFT);
		resHapmapLabel.setText("HapMap data set:");
		
		res_hapmap = new Text(ref_res_group,SWT.BORDER);
		res_hapmap.setText(RES_HAPMAP);
		res_hapmap.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		
		Button browseResHapmap = new Button(ref_res_group, SWT.RIGHT);
		browseResHapmap.setText("Browse");
		
		RES_OMNI = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.RES_OMNI);

		Label resOmniLabel = new Label(ref_res_group,SWT.LEFT);
		resOmniLabel.setText("Omni data set:");
		
		res_omni = new Text(ref_res_group,SWT.BORDER);
		res_omni.setText(RES_OMNI);
		res_omni.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		
		Button browseResOmni = new Button(ref_res_group, SWT.RIGHT);
		browseResOmni.setText("Browse");
		
		RES_1000G_SNPS = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.RES_1000G_SNPS);
		
		Label res1000GsnpsLabel = new Label(ref_res_group,SWT.LEFT);
		res1000GsnpsLabel.setText("1000G SNPs data set:");
		
		res_1000G_SNPS = new Text(ref_res_group,SWT.BORDER);
		res_1000G_SNPS.setText(RES_1000G_SNPS);
		res_1000G_SNPS.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		
		Button browseRes1000Gsnps = new Button(ref_res_group, SWT.RIGHT);
		browseRes1000Gsnps.setText("Browse");
		
		RES_1000G_INDELS = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.RES_1000G_INDELS);
		
		Label res1000GIndelsLabel = new Label(ref_res_group,SWT.LEFT);
		res1000GIndelsLabel.setText("1000G Indels data set:");
		
		res_1000G_Indels = new Text(ref_res_group,SWT.BORDER);
		res_1000G_Indels.setText(RES_1000G_INDELS);
		res_1000G_Indels.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		
		Button browseRes1000GIndels = new Button(ref_res_group, SWT.RIGHT);
		browseRes1000GIndels.setText("Browse");
		
		RES_DBSNP = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.RES_DBSNP);
		
		Label resDbsnpLabel = new Label(ref_res_group,SWT.LEFT);
		resDbsnpLabel.setText("dbSNP data set:");
		
		res_dbsnp = new Text(ref_res_group,SWT.BORDER);
		res_dbsnp.setText(RES_DBSNP);
		res_dbsnp.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		
		Button browseResDbsnp = new Button(ref_res_group, SWT.RIGHT);
		browseResDbsnp.setText("Browse");
		
		RES_MILLS = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.RES_MILLS);
		
		Label resMillsLabel = new Label(ref_res_group,SWT.LEFT);
		resMillsLabel.setText("Mills data set:");
		
		res_mills = new Text(ref_res_group,SWT.BORDER);
		res_mills.setText(RES_MILLS);
		res_mills.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		
		Button browseResMills = new Button(ref_res_group, SWT.RIGHT);
		browseResMills.setText("Browse");
		
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
		
		USE_HTE = IBISKNIMENodesPlugin.getBooleanPreference(IBISKNIMENodesPlugin.USE_HTE);
		checkHTE = new Button(use_hte,SWT.CHECK);
		checkHTE.setText("Use HTE?");
		checkHTE.setSelection(USE_HTE);
		
		
		Label thresholdLabel = new Label(use_hte,SWT.RIGHT);
		thresholdLabel.setText("Global threshold:");
		
		THRESHOLD = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.THRESHOLD);

		thresholdText = new Text(use_hte, SWT.BORDER);
		thresholdText.setText(THRESHOLD);
		thresholdText.setEnabled(USE_HTE);
		
		Label dbFileLabel = new Label(use_hte, SWT.LEFT);
		dbFileLabel.setText("Use existing db file:");
		
		DB_FILE = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.DB_FILE);
		
		dbFile = new Text(use_hte,SWT.BORDER);
		dbFile.setText(DB_FILE);
		dbFile.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		dbFile.setEnabled(USE_HTE);
		dbFile.setEditable(false);
		
		Button browseDBFile = new Button(use_hte, SWT.RIGHT);
		browseDBFile.setText("Browse");
		browseDBFile.setEnabled(USE_HTE);
		
		Button createDBFile = new Button(use_hte, SWT.NONE);
		createDBFile.setText("Create new db file");
		createDBFile.setEnabled(USE_HTE);
		
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
		
		NOTIFY = IBISKNIMENodesPlugin.getBooleanPreference(IBISKNIMENodesPlugin.NOTIFY);
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
		
		EMAILHOST = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.EMAIL_HOST);
		
		Label emailHost = new Label(email, SWT.LEFT);
		emailHost.setText("Email host:");
		email_host = new Text(email,SWT.BORDER);
		email_host.setText(EMAILHOST);
		email_host.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		email_host.setEnabled(NOTIFY);
		
		EMAILSENDER = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.EMAIL_SENDER);
		
		Label emailSender = new Label(email, SWT.LEFT);
		emailSender.setText("Email sender:");
		email_sender = new Text(email,SWT.BORDER);
		email_sender.setText(EMAILSENDER);
		email_sender.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		email_sender.setEnabled(NOTIFY);
		
		EMAILRECEIVER = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.EMAIL_RECEIVER);
		
		Label emailReceiver = new Label(email, SWT.LEFT);
		emailReceiver.setText("Email address (send to):");
		email_receiver = new Text(email,SWT.BORDER);
		email_receiver.setText(EMAILRECEIVER);
		email_receiver.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		email_receiver.setEnabled(NOTIFY);

		
		/* event handles */
		browseSearchDir.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				browseSearchDir.setEnabled(false);
				cancelSearch.setEnabled(true);
				String path = getDirPath(shell, "Select search directory");
				IBISKNIMENodesPlugin.getDefault().startSearchThread(path, table, browseSearchDir, cancelSearch);
			}
		});
		
		cancelSearch.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				IBISKNIMENodesPlugin.getDefault().cancelSearchThread();
				cancelSearch.setEnabled(false);
				browseSearchDir.setEnabled(true);
			}
		});
		
		downloader.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				String path = getDirPath(shell, "Select download directory");
				downloadBinaries(path);
			}
		});
		
		edit.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				
				if(table.getSelectionCount()==0) {
					return;
				}
				TableItem item = table.getSelection()[0];
				
				String file = getFilePath(shell, "Select directory path to "+item.getText(0),null, "");
				editBinary(item, file);
			}
		});
	
		
		browseRefGenome.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				String file = getFilePath(shell, "Select reference genome","*.fa;*.fasta", refGenome.getText());
				REF_GENOME = file;
				refGenome.setText(REF_GENOME);
			}
		});
		
		browseResHapmap.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				String file = getFilePath(shell, "Select HapMap dataset","*.vcf", res_hapmap.getText());
				RES_HAPMAP = file;
				res_hapmap.setText(RES_HAPMAP);
			}
		});
		
		browseResOmni.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				String file = getFilePath(shell, "Select Omni dataset","*.vcf", res_omni.getText());
				RES_OMNI = file;
				res_omni.setText(RES_OMNI);
			}
		});
		
		browseRes1000Gsnps.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				String file = getFilePath(shell, "Select 1000G SNPs data set","*.vcf", res_1000G_SNPS.getText());
				RES_1000G_SNPS = file;
				res_1000G_SNPS.setText(RES_1000G_SNPS);
			}
		});
		
		browseRes1000GIndels.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				String file = getFilePath(shell, "Select 1000G Indels data set","*.vcf", res_1000G_Indels.getText());
				RES_1000G_INDELS = file;
				res_1000G_Indels.setText(RES_1000G_INDELS);
			}
		});
		
		browseResDbsnp.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				String file = getFilePath(shell, "Select dbSNP dataset","*.vcf", res_dbsnp.getText());
				RES_DBSNP = file;
				res_dbsnp.setText(RES_DBSNP);
			}
		});
		
		browseResMills.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				String file = getFilePath(shell, "Select Mills dataset","*.vcf", res_mills.getText());
				RES_MILLS = file;
				res_mills.setText(RES_MILLS);
			}
		});
		
		browseDBFile.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				String file = getFilePath(shell, "Select HTE database file","*.db;*.sql", "");
				selectDBFile(file);
			}
		});
		
		createDBFile.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				String path = getDirPath(shell,"Choose directory in which database file will be stored");
				createDBFile(path);
			}
		});
		
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
		
		IBISKNIMENodesPlugin.setAllFieldsToDefault();
		
		for(TableItem i: table.getItems()) {
			i.setText(1,"");
		}
		
		checkOverwrite.setSelection(IBISKNIMENodesPlugin.OVERWRITE_DEFAULT);
		checkHTE.setSelection(IBISKNIMENodesPlugin.HTE_DEFAULT);
		thresholdText.setText(IBISKNIMENodesPlugin.THRESHOLD_DEFAULT+"");
		dbFile.setText("");
		checkNotify.setSelection(IBISKNIMENodesPlugin.NOTIFY_DEFAULT);
		refGenome.setText("");
		res_hapmap.setText("");
		res_omni.setText("");
		res_1000G_SNPS.setText("");
		res_1000G_Indels.setText("");
		res_dbsnp.setText("");
		res_mills.setText("");
		email_host.setText(IBISKNIMENodesPlugin.EMAIL_HOST_DEFAULT);
		email_sender.setText(IBISKNIMENodesPlugin.EMAIL_SENDER_DEFAULT);
		email_receiver.setText(IBISKNIMENodesPlugin.EMAIL_RECEIVER_DEFAULT);
	}
	
	/*
	 * The user has pressed Ok or Apply. Store/apply 
	 * this page's values appropriately.
	 */	
	public boolean performOk() {
		
		OVERWRITE = checkOverwrite.getSelection();
		IBISKNIMENodesPlugin.setBooleanPreference(IBISKNIMENodesPlugin.OVERWRITE, OVERWRITE);
		
		try {
				REF_GENOME = IO.processFilePath(refGenome.getText());
	
			if(!REF_GENOME.equals("") && !FileValidator.checkFastaFormat(REF_GENOME)) {
				JOptionPane.showMessageDialog(null,
						"Reference (genome) sequence file is not in FastA format or does not contain nucleotide sequences!",
					    "Error",
					    JOptionPane.ERROR_MESSAGE);
				return false;
			}
			IBISKNIMENodesPlugin.setStringPreference(IBISKNIMENodesPlugin.REF_GENOME, REF_GENOME);
			
			RES_HAPMAP = IO.processFilePath(res_hapmap.getText());
			IBISKNIMENodesPlugin.setStringPreference(IBISKNIMENodesPlugin.RES_HAPMAP, RES_HAPMAP);
			
			RES_OMNI = IO.processFilePath(res_omni.getText());
			IBISKNIMENodesPlugin.setStringPreference(IBISKNIMENodesPlugin.RES_OMNI, RES_OMNI);
			
			RES_1000G_SNPS = IO.processFilePath(res_1000G_SNPS.getText());
			IBISKNIMENodesPlugin.setStringPreference(IBISKNIMENodesPlugin.RES_1000G_SNPS, RES_1000G_SNPS);
			
			RES_1000G_INDELS = IO.processFilePath(res_1000G_Indels.getText());
			IBISKNIMENodesPlugin.setStringPreference(IBISKNIMENodesPlugin.RES_1000G_INDELS, RES_1000G_INDELS);
			
			RES_DBSNP = IO.processFilePath(res_dbsnp.getText());
			IBISKNIMENodesPlugin.setStringPreference(IBISKNIMENodesPlugin.RES_DBSNP, RES_DBSNP);
			
			RES_MILLS = IO.processFilePath(res_mills.getText());
			IBISKNIMENodesPlugin.setStringPreference(IBISKNIMENodesPlugin.RES_MILLS, RES_MILLS);
		
		} catch (InvalidSettingsException e1) {
		}
		
		
		USE_HTE = checkHTE.getSelection();
		IBISKNIMENodesPlugin.setBooleanPreference(IBISKNIMENodesPlugin.USE_HTE, USE_HTE);
		
		if(!USE_HTE) {
			IBISKNIMENodesPlugin.setBooleanPreference(IBISKNIMENodesPlugin.NOTIFY, USE_HTE);
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
		IBISKNIMENodesPlugin.setStringPreference(IBISKNIMENodesPlugin.THRESHOLD, THRESHOLD);
		
		try {
			DB_FILE = IO.processFilePath(dbFile.getText());
		} catch (InvalidSettingsException e1) {
		}
		
		if(DB_FILE.equals("")) {
			JOptionPane.showMessageDialog(null,
				    "HTE requires a SQLite database. Please create a new one or choose an existing one.",
				    "Error",
				    JOptionPane.ERROR_MESSAGE);
			return false;
		}
		IBISKNIMENodesPlugin.setStringPreference(IBISKNIMENodesPlugin.DB_FILE, DB_FILE);		

		EMAILHOST = email_host.getText();
		IBISKNIMENodesPlugin.setStringPreference(IBISKNIMENodesPlugin.EMAIL_HOST, EMAILHOST);
		
		EMAILSENDER = email_sender.getText();
		IBISKNIMENodesPlugin.setStringPreference(IBISKNIMENodesPlugin.EMAIL_SENDER, EMAILSENDER);
		
		EMAILRECEIVER = email_receiver.getText();
		IBISKNIMENodesPlugin.setStringPreference(IBISKNIMENodesPlugin.EMAIL_RECEIVER, EMAILRECEIVER);
		
		NOTIFY = checkNotify.getSelection();
		IBISKNIMENodesPlugin.setBooleanPreference(IBISKNIMENodesPlugin.NOTIFY, NOTIFY);
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
	
	private void downloadBinaries(String dir){
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
			
		for (String tool : IBISKNIMENodesPlugin.TOOLS.keySet()) {
			if (IBISKNIMENodesPlugin.TOOLS.get(tool)) {
				if (IBISKNIMENodesPlugin.getStringPreference(tool).equals("")) {
					try {
						File f = new File(dir + File.separatorChar + tool);
						if (!f.exists()) {
							FileUtils.copyURLToFile(new URL(DOWNLOAD_PATH + tool), f);
							f.setExecutable(true, false);
							IBISKNIMENodesPlugin.setStringPreference(tool, dir + "/"+ tool);
							
							HashSet<String> deps = DEPENDENCIES.get(tool);
							if(deps == null) continue;
							for(String dep: deps) {
								File f1 = new File(dir + File.separatorChar + dep);
								if (!f1.exists()) {
									FileUtils.copyURLToFile(new URL(DOWNLOAD_PATH + dep), f1);
									f1.setExecutable(true, false);
								}
							}
						}
					} catch (IOException e) {
						LOGGER.error("Downloading " + tool+ " failed! Message: " + e.getMessage());
					}
				}
			}
		}
		
		for(TableItem i: table.getItems()) {
			i.setText(1,IBISKNIMENodesPlugin.getStringPreference(i.getText(0)));
		}
	}
	

	
	private String getFilePath(Shell shell, String text, String ext, String content) {
		FileDialog fdl = new FileDialog(shell);
		fdl.setText(text);
		if(ext == null) {
			fdl.setFilterExtensions(new String[]{ext});
		}
		if(content.equals("")) {
			fdl.setFilterPath(System.getProperty("user.home"));
		} else {
			fdl.setFilterPath(new File(content).getParent());
		}
		
		String path = fdl.open();
		return path;
	}
	
	private String getDirPath(Shell shell, String text) {
		DirectoryDialog dl = new DirectoryDialog(shell);
		dl.setText(text);
		dl.setFilterPath("~/");
		String path = dl.open();
		return path;
	}
	
	private void editBinary(TableItem item, String file) {
		if(file==null) return;
		
		try {
			CheckUtils.checkSourceFile(file);
			IBISKNIMENodesPlugin.setStringPreference(item.getText(0), file);
			item.setText(1,file);
		} catch (InvalidSettingsException e) {
			JOptionPane.showMessageDialog(null,
				    "This file cannot be used: "+e.getMessage(),
				    "Error",
				    JOptionPane.ERROR_MESSAGE);
			return;
		}
	}
	
	private void selectDBFile(String path) {
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
	
	public void createDBFile(String dir) {
		if(dir==null) return;
		try {
			if(Files.exists(Paths.get(dir+"/hte.db"))) {
				String e = "Database hte.db exists in "+dir+". Please use another location or rename the existing database.";
				throw new SQLException(e);
			}
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

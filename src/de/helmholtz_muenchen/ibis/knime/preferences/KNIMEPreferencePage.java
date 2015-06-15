
package de.helmholtz_muenchen.ibis.knime.preferences;


import java.io.File;
import java.io.IOException;
import java.net.URL;

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
import org.eclipse.swt.widgets.Group;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Shell;
import org.eclipse.swt.widgets.Text;
import org.eclipse.ui.IWorkbench;
import org.eclipse.ui.IWorkbenchPreferencePage;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;

/**
 * @author hastreiter
 */
public class KNIMEPreferencePage extends PreferencePage implements
        IWorkbenchPreferencePage {

	public static String TOOL_LOCATION;
	public static boolean USE_HTE;
	public static String THRESHOLD;
	
	private Text binsDirectory;
	private Text thresholdText;


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
		

		TOOL_LOCATION = IBISKNIMENodesPlugin.getDefault().getToolDirPreference();
		
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
		
		binsDirectory = new Text(downloadBins, SWT.BORDER);
		binsDirectory.setText(TOOL_LOCATION);
		
		
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
		
		Composite use_hte = new Composite(htePrefs,SWT.NONE);
		use_hte.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		GridLayout hteLayout = new GridLayout();
		hteLayout.numColumns = 1;
		use_hte.setLayout(hteLayout);
		
		USE_HTE = IBISKNIMENodesPlugin.getDefault().getHTEPreference();
		Button checkHTE = new Button(use_hte,SWT.CHECK);
		checkHTE.setText("Use HTE?");
		checkHTE.setSelection(USE_HTE);
		checkHTE.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				Button btn = (Button) e.getSource();
				USE_HTE = btn.getSelection();
			}
		});
		
		Composite dbPrefs = new Composite(htePrefs,SWT.NONE);
		dbPrefs.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		GridLayout dbLayout = new GridLayout();
		dbLayout.numColumns = 2;
		dbPrefs.setLayout(dbLayout);
		
		Label thresholdLabel = new Label(dbPrefs,SWT.NONE);
		thresholdLabel.setText("Global threshold:");
		
		THRESHOLD = IBISKNIMENodesPlugin.getDefault().getThresholdPreference();

		thresholdText = new Text(dbPrefs, SWT.BORDER);
		thresholdText.setText(THRESHOLD);
		
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
		
		THRESHOLD = thresholdText.getText();
		iknp.setThresholdPreference(THRESHOLD);
		System.out.println("Setting THREHOLD to: "+THRESHOLD);
		
		return super.performOk();
	}
	
	private void downloadBinaries(){
		
		String dir = TOOL_LOCATION;
		System.out.println(dir);

		if(dir!=null){
		try {
			File bwa = new File(dir+"/bwa");
			FileUtils.copyURLToFile(new URL("ftp://ftpmips.helmholtz-muenchen.de/Incoming/KNIME_BIN/bwa"), bwa );
			bwa.setExecutable(true,false);
			
			File pindel = new File(dir+"/pindel");
			FileUtils.copyURLToFile(new URL("ftp://ftpmips.helmholtz-muenchen.de/Incoming/KNIME_BIN/pindel"), pindel);
			pindel.setExecutable(true,false);
			
			File pindel2vcf = new File(dir+"/pindel2vcf");
			FileUtils.copyURLToFile(new URL("ftp://ftpmips.helmholtz-muenchen.de/Incoming/KNIME_BIN/pindel2vcf"), pindel2vcf);
			pindel2vcf.setExecutable(true,false);
			
			FileUtils.copyURLToFile(new URL("ftp://ftpmips.helmholtz-muenchen.de/Incoming/KNIME_BIN/GenomeAnalysisTK.jar"), new File(dir+"/GenomeAnalysisTK.jar"));

		} catch (IOException e) {
			e.printStackTrace();
		}
		}
	}
	
	
	private void selectToolFolder(Shell shell){
		
		DirectoryDialog dlg = new DirectoryDialog(shell);
		dlg.setText("Choose directory in which tool binaries will be stored");
		dlg.setFilterPath("~/");
		String dir = dlg.open();
		System.out.println(dir);
		TOOL_LOCATION=dir;
		binsDirectory.setText(TOOL_LOCATION);
	}
	
}

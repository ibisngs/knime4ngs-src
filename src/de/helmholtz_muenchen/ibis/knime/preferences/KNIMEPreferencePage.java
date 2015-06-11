
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




	public static String TOOL_LOCATION = null;
	private Text textField;
	private Text textField2;


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

		// Sets the layout data for the top composite's 
		// place in its parent's layout.
		top.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));

		// Sets the layout for the top composite's 
		// children to populate.
		top.setLayout(new GridLayout());
				
		
		TOOL_LOCATION = IBISKNIMENodesPlugin.getDefault().getToolDirPreference();
		
		
		Composite DownloadGroup = new Composite(top, SWT.NONE);
		DownloadGroup.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		GridLayout DownloadLayout = new GridLayout();
		DownloadLayout.numColumns = 2;
		DownloadLayout.marginHeight = 0;
		DownloadLayout.marginWidth = 0;
		DownloadGroup.setLayout(DownloadLayout);
		
		// Create a composite for the add and remove buttons.
		Composite buttonGroup = new Composite(DownloadGroup, SWT.NONE);
		buttonGroup.setLayoutData(new GridData());
		GridLayout buttonLayout = new GridLayout();
		buttonLayout.marginHeight = 0;
		buttonLayout.marginWidth = 0;
		buttonGroup.setLayout(buttonLayout);

		Button addTag = new Button(buttonGroup, SWT.NONE);
		addTag.setText("Download Tool Binaries");
		final Shell shell = new Shell(parent.getDisplay());
		addTag.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				downloadBinaries(shell);
			}
		});
		addTag.setSize(200, 20);
		
		textField = new Text(DownloadGroup, SWT.BORDER);
		
		GridData textData = new GridData(GridData.FILL_HORIZONTAL);
		textData.verticalAlignment = GridData.BEGINNING;
		textField.setLayoutData(textData);
		textField.setText(TOOL_LOCATION);
		

		
		
		Composite SelectExistingGroup = new Composite(top, SWT.NONE);
		SelectExistingGroup.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		GridLayout SelectExistingLayout = new GridLayout();
		SelectExistingLayout.numColumns = 2;
		SelectExistingLayout.marginHeight = 0;
		SelectExistingLayout.marginWidth = 0;
		SelectExistingGroup.setLayout(SelectExistingLayout);
		
		// Create a composite for the add and remove buttons.
		Composite buttonGroup2 = new Composite(SelectExistingGroup, SWT.NONE);
		buttonGroup2.setLayoutData(new GridData());
		GridLayout buttonLayout2 = new GridLayout();
		buttonLayout2.marginHeight = 0;
		buttonLayout2.marginWidth = 0;
		buttonGroup2.setLayout(buttonLayout2);

		Button addTag2 = new Button(buttonGroup2, SWT.NONE);
		addTag2.setText("Or select existing folder with binaries");
		final Shell shell2 = new Shell(parent.getDisplay());
		addTag2.addSelectionListener(new SelectionAdapter() {
			public void widgetSelected(SelectionEvent e) {
				selectToolFolder(shell2);
			}
		});
		addTag2.setSize(200, 20);
		
		textField2 = new Text(SelectExistingGroup, SWT.BORDER);
		
		GridData textData2 = new GridData(GridData.FILL_HORIZONTAL);
		textData2.verticalAlignment = GridData.BEGINNING;
		textField2.setLayoutData(textData2);
		textField2.setText(TOOL_LOCATION);
		


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
		IBISKNIMENodesPlugin.getDefault().
		setToolDirPreference(TOOL_LOCATION);
		System.out.println("Seeting TOOL_LOCATION to: "+TOOL_LOCATION);
		return super.performOk();
	}
	
	
	private void downloadBinaries(Shell shell){
		
		DirectoryDialog dlg = new DirectoryDialog(shell);
		dlg.setText("Choose directory in which tool binaries will be stored");
		dlg.setFilterPath("~/");
		String dir = dlg.open();
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

			TOOL_LOCATION=dir;
			textField.setText(TOOL_LOCATION);
		} catch (IOException e) {
			// TODO Auto-generated catch block
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
		textField.setText(TOOL_LOCATION);
		textField2.setText(TOOL_LOCATION);
	}
	
}

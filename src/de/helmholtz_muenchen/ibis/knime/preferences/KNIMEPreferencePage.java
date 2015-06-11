
package de.helmholtz_muenchen.ibis.knime.preferences;


import java.io.File;
import java.io.IOException;
import java.net.URL;

import org.apache.commons.io.FileUtils;
//import org.eclipse.jface.dialogs.IDialogConstants;
import org.eclipse.jface.preference.FieldEditor;
import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.jface.preference.PreferencePage;
import org.eclipse.jface.preference.RadioGroupFieldEditor;
import org.eclipse.swt.SWT;
import org.eclipse.swt.events.SelectionAdapter;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.layout.GridData;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Control;
import org.eclipse.swt.widgets.DirectoryDialog;
//import org.eclipse.swt.widgets.Display;
//import org.eclipse.swt.widgets.Event;
//import org.eclipse.swt.widgets.Label;
//import org.eclipse.swt.widgets.List;
//import org.eclipse.swt.widgets.Listener;
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
		
		
		Composite addRemoveGroup = new Composite(top, SWT.NONE);
		addRemoveGroup.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		GridLayout addRemoveLayout = new GridLayout();
		addRemoveLayout.numColumns = 2;
		addRemoveLayout.marginHeight = 0;
		addRemoveLayout.marginWidth = 0;
		addRemoveGroup.setLayout(addRemoveLayout);
		
		// Create a composite for the add and remove buttons.
		Composite buttonGroup = new Composite(addRemoveGroup, SWT.NONE);
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
		
		textField = new Text(addRemoveGroup, SWT.BORDER);
		
		GridData textData = new GridData(GridData.FILL_HORIZONTAL);
		textData.verticalAlignment = GridData.BEGINNING;
		textField.setLayoutData(textData);
		textField.setText(TOOL_LOCATION);
		

		
		
		
		
		Composite addRemoveGroup2 = new Composite(top, SWT.NONE);
		addRemoveGroup2.setLayoutData(new GridData(GridData.FILL_HORIZONTAL));
		GridLayout addRemoveLayout2 = new GridLayout();
		addRemoveLayout2.numColumns = 2;
		addRemoveLayout2.marginHeight = 0;
		addRemoveLayout2.marginWidth = 0;
		addRemoveGroup2.setLayout(addRemoveLayout2);
		
		// Create a composite for the add and remove buttons.
		Composite buttonGroup2 = new Composite(addRemoveGroup2, SWT.NONE);
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
		
		textField2 = new Text(addRemoveGroup2, SWT.BORDER);
		
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
//		errors.loadDefault();
//		exemptTagsList.setItems(
//				IBISKNIMENodesPlugin.getDefault().
//				getDefaultExemptTagsPreference());
		// getDefaultExemptTagsPreference() is a convenience
		// method which retrieves the default preference from
		// the preference store.
		super.performDefaults();
	}
	
	/*
	 * The user has pressed Ok or Apply. Store/apply 
	 * this page's values appropriately.
	 */	
	public boolean performOk() {
//		errors.store();
//		IBISKNIMENodesPlugin.getDefault().
//			setExemptTagsPreference(exemptTagsList.getItems());
		IBISKNIMENodesPlugin.getDefault().
		setToolDirPreference(TOOL_LOCATION);
		System.out.println("Seeting TOOL_LOCATION to: "+TOOL_LOCATION);
		return super.performOk();
	}
	
//	private void addTag() {
//		String tag = textField.getText();
//		if (tag != null && tag.length() > 0)
//			exemptTagsList.add(tag);
//		textField.setText("");
//	}
	
	/*
	 * Sets the enablement of the remove button depending
	 * on the selection in the list.
	 */
//	private void selectionChanged() {
//		int index = exemptTagsList.getSelectionIndex();
//		removeTag.setEnabled(index >= 0);		
//	}
	
	
	private void downloadBinaries(Shell shell){
		
		DirectoryDialog dlg = new DirectoryDialog(shell);
		dlg.setText("Choose directory in which tool binaries will be stored");
		dlg.setFilterPath("~/");
		String dir = dlg.open();
		System.out.println(dir);

		if(dir!=null){
		try {
//			URL url = new URL("ftp://ftpmips.helmholtz-muenchen.de/Incoming/KNIME_BIN/");
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
	}
	
}

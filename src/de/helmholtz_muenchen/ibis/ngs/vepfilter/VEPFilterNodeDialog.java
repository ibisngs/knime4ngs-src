package de.helmholtz_muenchen.ibis.ngs.vepfilter;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentButton;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringListSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;


/**
 * <code>NodeDialog</code> for the "VEPFilter" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author tim.jeske
 */
public class VEPFilterNodeDialog extends HTExecutorNodeDialog {
	
	private DialogComponentStringListSelection DC_TERM_DISPLAY;
	private SettingsModelBoolean my_overwrite;
    
    public void addToolDialogComponents() {
    	
    	final SettingsModelString vep_script = new SettingsModelString(VEPFilterNodeModel.CFGKEY_VEP_SCRIPT,"");
    	final SettingsModelOptionalString filter = new SettingsModelOptionalString(VEPFilterNodeModel.CFGKEY_FILTER,"",false);
    	final SettingsModelString outfolder = new SettingsModelString(VEPFilterNodeModel.CFGKEY_OUTFOLDER,"");
    	final SettingsModelBoolean overwrite = new SettingsModelBoolean(VEPFilterNodeModel.CFGKEY_OVERWRITE,false);
    	my_overwrite = overwrite;
    	
    	final SettingsModelStringArray chosen_terms = new SettingsModelStringArray(VEPFilterNodeModel.CFGKEY_TERM_LIST, new String[]{VEPFilterNodeModel.DEFAULT_TERM});
    	
    	DC_TERM_DISPLAY = new DialogComponentStringListSelection(chosen_terms, "", VEPFilterNodeModel.SO_TERMS);

    	DialogComponentButton CHOOSE_LOF_TERMS = new DialogComponentButton("Choose LOF terms");
    	DialogComponentOptionalString DC_FILTER = new DialogComponentOptionalString(filter,"Conditions");
    	
    	addPrefPageSetting(vep_script, IBISKNIMENodesPlugin.VEP_FILTER);
    	
    	createNewGroup("Filter by Consequence");
		DC_TERM_DISPLAY.setVisibleRowCount(20);
		addDialogComponent(DC_TERM_DISPLAY);
		
		CHOOSE_LOF_TERMS.setToolTipText("Default SO terms comprise variants that are considered as LoF variants.");
		addDialogComponent(CHOOSE_LOF_TERMS);
		
		createNewGroup("Further filtering conditions");
		DC_FILTER.setToolTipText("Define additional comma separated filters,e.g.\"LoF is HC,ExAC_AF>0.05");
		addDialogComponent(DC_FILTER);

		
		CHOOSE_LOF_TERMS.addActionListener(new ActionListener () {
			@Override
			public void actionPerformed(ActionEvent e) {
				chosen_terms.setStringArrayValue(VEPFilterNodeModel.LOF_TERMS);
			}
		});
	
		createNewGroup("Folder for output files");
    	addDialogComponent(new DialogComponentFileChooser(outfolder, "his_id_VEP_Filtered_OUT", 0, true));
    	addDialogComponent(new DialogComponentBoolean(overwrite, "Overwrite, if output files exist?"));
    }
    
    public void onOpen() {
    	super.onOpen();
    	boolean use_hte = IBISKNIMENodesPlugin.getBooleanPreference(IBISKNIMENodesPlugin.USE_HTE);
    	
    	if(use_hte) {
    		my_overwrite.setBooleanValue(true);
    		my_overwrite.setEnabled(false);
    	} else {
//    		my_overwrite.setBooleanValue(false);
    		my_overwrite.setEnabled(true);
    	}
    }
}
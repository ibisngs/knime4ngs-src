package de.helmholtz_muenchen.ibis.ngs.fastqcStatisticMerger;

import java.util.HashMap;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * Dialog for FastQC Statistic Merger Node 
 * @author Michael Kluge
 *
 */
public class FastqcStatisticMergerNodeDialog extends DefaultNodeSettingsPane {

    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_INPUT_FOLDER	= new SettingsModelString(FastqcStatisticMergerNodeModel.CFGKEY_INPUT_FOLDER, FastqcStatisticMergerNodeModel.DEFAULT_INPUT_FOLDER);
    private final SettingsModelString SET_OUTPUT_FOLDER = new SettingsModelString(FastqcStatisticMergerNodeModel.CFGKEY_OUTPUT_FOLDER, FastqcStatisticMergerNodeModel.DEFAULT_OUTPUT_FOLDER);
    private final SettingsModelBoolean SET_ALL_MODULES = new SettingsModelBoolean(FastqcStatisticMergerNodeModel.CFGKEY_MODULE_ALL, false);
    
    // store boolean settings
    private final HashMap<String, SettingsModelBoolean> MODULE_OPTIONS = new HashMap<String, SettingsModelBoolean>();
    
    protected FastqcStatisticMergerNodeDialog() {
    	// create open file/folder components
        DialogComponentFileChooser dcInputFolder 	= new DialogComponentFileChooser(SET_INPUT_FOLDER, "his_id_INPUT_FOLDER" + this.getClass().getCanonicalName(), 0, true);
        DialogComponentFileChooser dcOutputFolder 	= new DialogComponentFileChooser(SET_OUTPUT_FOLDER, "his_id_OUTPUT_FOLDER" + this.getClass().getCanonicalName(), 0, true);
        
       	// set a new title to them+
       	dcInputFolder.setBorderTitle("path to folder with FastQC Outputs");
       	dcOutputFolder.setBorderTitle("path to output folder");
     
       	// add groups and components
        createNewGroup("input");
        addDialogComponent(dcInputFolder);
        
        createNewGroup("output");
        addDialogComponent(dcOutputFolder);
        
        // add the module names
        createNewGroup("modules to merge");
        
        for(String moduleName : FastqcStatisticMergerNodeModel.MODULE_NAMES) {
        	SettingsModelBoolean option = new SettingsModelBoolean(moduleName, false);
        	DialogComponentBoolean com = new DialogComponentBoolean(option, moduleName);
        	addDialogComponent(com);
        	MODULE_OPTIONS.put(moduleName, option);
        }
   
        // add all modules
        DialogComponentBoolean com = new DialogComponentBoolean(SET_ALL_MODULES, FastqcStatisticMergerNodeModel.CFGKEY_MODULE_ALL);
    	addDialogComponent(com);
        
        // add change listener to all checkbox
        SET_ALL_MODULES.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				for(String moduleName : FastqcStatisticMergerNodeModel.MODULE_NAMES) {
					SettingsModelBoolean option = MODULE_OPTIONS.get(moduleName);

					if(SET_ALL_MODULES.getBooleanValue())
						option.setBooleanValue(true); // enable checkbox
					else 
						option.setBooleanValue(false); // disable checkbox
				}
			}
        });
    }
}


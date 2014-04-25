package de.helmholtz_muenchen.ibis.ngs.fastqcStatisticMerger;

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
    private final SettingsModelString SET_INPUT_FOLDER	= FastqcStatisticMergerNodeModel.getSettingsModelString(FastqcStatisticMergerNodeModel.CFGKEY_INPUT_FOLDER, this);
    private final SettingsModelString SET_OUTPUT_FOLDER = FastqcStatisticMergerNodeModel.getSettingsModelString(FastqcStatisticMergerNodeModel.CFGKEY_OUTPUT_FOLDER, this);
    private final SettingsModelBoolean SET_ALL_MODULES = FastqcStatisticMergerNodeModel.getSettingsModelBoolean(FastqcStatisticMergerNodeModel.CFGKEY_MODULE_ALL, this);
    private final FastqcStatisticMergerNodeDialog THIS;
    
    protected FastqcStatisticMergerNodeDialog() {
    	THIS = this; // is used in the eventlistener because of the hascode
    	
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
        	SettingsModelBoolean option = FastqcStatisticMergerNodeModel.getSettingsModelBoolean(moduleName, this);
        	DialogComponentBoolean com = new DialogComponentBoolean(option, moduleName);
        	addDialogComponent(com);
        }
        
        // add change listener to all checkbox
        SET_ALL_MODULES.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				for(String moduleName : FastqcStatisticMergerNodeModel.MODULE_NAMES) {
					SettingsModelBoolean option = FastqcStatisticMergerNodeModel.getSettingsModelBoolean(moduleName, THIS);
					
					if(SET_ALL_MODULES.getBooleanValue())
						option.setBooleanValue(true); // enable checkbox
					else 
						option.setBooleanValue(false); // disable checkbox
				}
			}
        });
    }
}


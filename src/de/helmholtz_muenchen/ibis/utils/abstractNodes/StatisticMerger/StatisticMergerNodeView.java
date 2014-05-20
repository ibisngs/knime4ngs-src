package de.helmholtz_muenchen.ibis.utils.abstractNodes.StatisticMerger;

import java.util.HashMap;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * Standard dialog for statistic mergers
 * @author Michael Kluge
 *
 */
public abstract class StatisticMergerNodeView extends DefaultNodeSettingsPane {

    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_INPUT_FOLDER	= new SettingsModelString(StatisticMergerNodeModel.CFGKEY_INPUT_FOLDER, StatisticMergerNodeModel.DEFAULT_INPUT_FOLDER);
    private final SettingsModelString SET_OUTPUT_FOLDER = new SettingsModelString(StatisticMergerNodeModel.CFGKEY_OUTPUT_FOLDER, StatisticMergerNodeModel.DEFAULT_OUTPUT_FOLDER);
    private final SettingsModelBoolean SET_ALL_MODULES = new SettingsModelBoolean(StatisticMergerNodeModel.CFGKEY_MODULE_ALL, false);
    
    // store boolean settings
    private final HashMap<String, SettingsModelBoolean> MODULE_OPTIONS = new HashMap<String, SettingsModelBoolean>();
    private final String MERGER_NAME;
    
    protected StatisticMergerNodeView() {
    	this.MERGER_NAME = this.getMergerName();
    	
    	// create open file/folder components
        DialogComponentFileChooser dcInputFolder 	= new DialogComponentFileChooser(SET_INPUT_FOLDER, "his_id_INPUT_FOLDER" + this.getClass().getCanonicalName(), 0, true);
        DialogComponentFileChooser dcOutputFolder 	= new DialogComponentFileChooser(SET_OUTPUT_FOLDER, "his_id_OUTPUT_FOLDER" + this.getClass().getCanonicalName(), 0, true);
        
       	// set a new title to them+
       	dcInputFolder.setBorderTitle("path to folder with " + this.MERGER_NAME + " Outputs");
       	dcOutputFolder.setBorderTitle("path to output folder");
     
       	// add groups and components
        createNewGroup("input");
        addDialogComponent(dcInputFolder);
        
        createNewGroup("output");
        addDialogComponent(dcOutputFolder);
        
        // add the module names
        createNewGroup("modules to merge");
        // add all modules
        for(String moduleName : StatisticMergerNodeModel.getModuleNames(this.MERGER_NAME)) {
        	SettingsModelBoolean option = new SettingsModelBoolean(moduleName, false);
        	DialogComponentBoolean com = new DialogComponentBoolean(option, moduleName);
        	addDialogComponent(com);
        	MODULE_OPTIONS.put(moduleName, option);
        }
   
        
    	// only, add this checkbox, if more than one element is in the list.
        if(StatisticMergerNodeModel.getModuleNames(this.MERGER_NAME).size() > 1) {
        	DialogComponentBoolean com = new DialogComponentBoolean(SET_ALL_MODULES, StatisticMergerNodeModel.CFGKEY_MODULE_ALL);
    		addDialogComponent(com);
        
	        // add change listener to all checkbox
	        SET_ALL_MODULES.addChangeListener(new ChangeListener() {
				@Override
				public void stateChanged(ChangeEvent arg0) {
					for(String moduleName : StatisticMergerNodeModel.getModuleNames(MERGER_NAME)) {
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
    
    /**
     * Name of the statistic merger
     * @return
     */
    public abstract String getMergerName();
}


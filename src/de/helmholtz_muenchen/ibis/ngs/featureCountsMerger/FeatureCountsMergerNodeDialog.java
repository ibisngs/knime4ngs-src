package de.helmholtz_muenchen.ibis.ngs.featureCountsMerger;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "FeatureCountsMerger" Node.
 * 
 * @author Michael Kluge
 */
public class FeatureCountsMergerNodeDialog extends DefaultNodeSettingsPane {

	private final SettingsModelString SET_OUTPUT_FILE = new SettingsModelString(FeatureCountsMergerNodeModel.CFGKEY_OUTPUT_FILE, FeatureCountsMergerNodeModel.DEFAULT_OUTPUT_FILE);
	private final SettingsModelBoolean SET_REMOVE_ENDING = new SettingsModelBoolean(FeatureCountsMergerNodeModel.CFGKEY_REMOVE_ENDING, FeatureCountsMergerNodeModel.DEFAULT_REMOVE_ENDING);
	private final SettingsModelBoolean SET_REMOVE_PATH = new SettingsModelBoolean(FeatureCountsMergerNodeModel.CFGKEY_REMOVE_PATH, FeatureCountsMergerNodeModel.DEFAULT_REMOVE_PATH);

    /**
     * New pane for configuring the FeatureCountsMerger node.
     */
    protected FeatureCountsMergerNodeDialog() {
		// create open file components
		DialogComponentFileChooser dcOutputFile = new DialogComponentFileChooser(SET_OUTPUT_FILE, "his_id_fcm_OUTPUT_FILE", 0, false);
		DialogComponentBoolean removePath 		= new DialogComponentBoolean(SET_REMOVE_PATH, "Remove the path of the file in the header");
		DialogComponentBoolean removeEnding 	= new DialogComponentBoolean(SET_REMOVE_ENDING, "Remove the ending of the file in the header");
		
		// set a new title to them
		dcOutputFile.setBorderTitle("Path to output file");
		 
		// add groups and components        
		createNewGroup("Output");
		addDialogComponent(dcOutputFile);
		
		createNewGroup("Header format options");
		addDialogComponent(removePath);
		addDialogComponent(removeEnding);
    }
}


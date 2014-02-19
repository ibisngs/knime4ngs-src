package de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentMultiLineString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>BinaryWrapperNodeDialog</code> provides some options to set a path to a binary and set additional commands.
 * @author Michael Kluge
 *
 */
public abstract class BinaryWrapperNodeDialog extends DefaultNodeSettingsPane {

    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_BINARY_PATH			= BinaryWrapperNodeModel.getSettingsModelString(BinaryWrapperNodeModel.CFGKEY_BINARY_PATH, this);
    private final SettingsModelString SET_ADDITIONAL_PARAMETER 	= BinaryWrapperNodeModel.getSettingsModelString(BinaryWrapperNodeModel.CFGKEY_ADDITIONAL_PARAMETER, this);
    private final SettingsModelString SET_PARAMETER_FILE		= BinaryWrapperNodeModel.getSettingsModelString(BinaryWrapperNodeModel.CFGKEY_PARAMETER_FILE, this);
    
	protected BinaryWrapperNodeDialog() {
		super();
		// rename default tab
		setDefaultTabTitle("general options (level 0 and 1)");
		
        // create open file/folder components
        DialogComponentFileChooser dcBinaryPath = new DialogComponentFileChooser(SET_BINARY_PATH, "his_id_BINARY_PATH", 0);
        DialogComponentFileChooser dcParameterFile 	= new DialogComponentFileChooser(SET_PARAMETER_FILE, "his_id_PARAMETER_FILE", 0);
        DialogComponentMultiLineString dcAdditionalParameter = new DialogComponentMultiLineString(SET_ADDITIONAL_PARAMETER, "additional parameters (level 1):", false, 45, 3);
        
        // set title
    	dcBinaryPath.setBorderTitle("path to " + getNameOfBinary() + " binary");
        dcParameterFile.setBorderTitle("path to parameter file (level 0)");
        
       	// add binary group and components
        createNewGroup("binary file path and additional parameters");
        addDialogComponent(dcBinaryPath);
        addDialogComponent(dcParameterFile);
        addDialogComponent(dcAdditionalParameter);
        
        // create new tab for GUI options
        createNewTab("level 2 options");
	}
	
	/**
	 * Should return the name of the binary for the GUI
	 * @return
	 */
	protected abstract String getNameOfBinary();
}

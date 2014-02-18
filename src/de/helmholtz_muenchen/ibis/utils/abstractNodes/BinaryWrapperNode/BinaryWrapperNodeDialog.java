package de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.ngs.star.StarNodeModel;

public abstract class BinaryWrapperNodeDialog extends DefaultNodeSettingsPane {

	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_BINARY_PATH = StarNodeModel.getSettingsModelString(BinaryWrapperNodeModel.CFGKEY_BINARY_PATH, this);
	
	protected BinaryWrapperNodeDialog() {
		super();
		
        // create open file/folder components
        DialogComponentFileChooser dcBinaryPath = new DialogComponentFileChooser(SET_BINARY_PATH, "his_id_BINARY_PATH", 0);
    	dcBinaryPath.setBorderTitle("path to STAR binary");
        
       	// add binary group and components
        createNewGroup("binary file and additional parameters");
        addDialogComponent(dcBinaryPath);
	}
}

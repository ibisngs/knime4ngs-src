package de.helmholtz_muenchen.ibis.ngs.fillmissinggenotypes;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "FillMissingGenotypes" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class FillMissingGenotypesNodeDialog extends DefaultNodeSettingsPane {

	private final SettingsModelString COVERAGEFOLDER 		= new SettingsModelString(FillMissingGenotypesNodeModel.CFGKEY_COVERAGEFILEFOLDER, "");
    private final SettingsModelString FILESUFFIX 			= new SettingsModelString(FillMissingGenotypesNodeModel.CFGKEY_FILESUFFIX, "");

	
	
    /**
     * New pane for configuring the FillMissingGenotypes node.
     */
    protected FillMissingGenotypesNodeDialog() {

    	createNewGroup("Folder in which Coverage files can be found");
    	addDialogComponent(new DialogComponentFileChooser(COVERAGEFOLDER, "Folder in which Coverage files can be found", 0,true,""));
    	createNewGroup("");
    	addDialogComponent(new DialogComponentString(FILESUFFIX, "File Suffix of Coverage Files which are used to fill ./. GTs"));
    	
    }
}


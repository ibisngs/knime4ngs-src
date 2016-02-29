package de.helmholtz_muenchen.ibis.ngs.gatkcombinegvcfs;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

/**
 * <code>NodeDialog</code> for the "CombineGVCFs" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class CombineGVCFsNodeDialog extends GATKNodeDialog {
	
	private SettingsModelString OUTFOLDER;

	@Override
	protected void addDialogComponent() {
		OUTFOLDER = new SettingsModelString(CombineGVCFsNodeModel.CFGKEY_OUTFOLDER, "");
		
		
	createNewGroup("Folder for output files");
    addDialogComponent(new DialogComponentFileChooser(OUTFOLDER, "his_id_VEP_OUT", 0, true));	
	}
}


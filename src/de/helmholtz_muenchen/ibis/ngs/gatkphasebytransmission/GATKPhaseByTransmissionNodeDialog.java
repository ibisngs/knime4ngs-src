package de.helmholtz_muenchen.ibis.ngs.gatkphasebytransmission;


import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

/**
 * <code>NodeDialog</code> for the "GATKPhaseByTransmission" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class GATKPhaseByTransmissionNodeDialog extends GATKNodeDialog {

	@Override
	protected void addDialogComponent() {
		SettingsModelString PED_FILE = new SettingsModelString(GATKPhaseByTransmissionNodeModel.CFGKEY_PED_FILE, "");
	    SettingsModelString DENOVO_PRIOR = new SettingsModelString(GATKPhaseByTransmissionNodeModel.CFGKEY_DENOVOPRIOR, "1.0E-8");
		
	    createNewTab("PhaseByTransmission");
	    createNewGroup("PED File");
    	DialogComponentFileChooser ped_file= new DialogComponentFileChooser(PED_FILE, "ped_file", 0, ".ped");
    	addDialogComponent(ped_file);
    	
    	createNewGroup("");
    	addDialogComponent(new DialogComponentString(DENOVO_PRIOR, "De Novo Prior"));
	}
}


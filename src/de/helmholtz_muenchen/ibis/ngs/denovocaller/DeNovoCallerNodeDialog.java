package de.helmholtz_muenchen.ibis.ngs.denovocaller;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "DeNovoCaller" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class DeNovoCallerNodeDialog extends DefaultNodeSettingsPane {

	/**
	 * Node Models
	 */
	private final SettingsModelString inputtype = new SettingsModelString(DeNovoCallerNodeModel.CFGKEY_INPUTTYPE,"");
	private final SettingsModelString pedfile = new SettingsModelString(DeNovoCallerNodeModel.CFGKEY_PEDFILE,"");
	
	
    /**
     * New pane for configuring the DeNovoCaller node.
     */
    protected DeNovoCallerNodeDialog() {

    	addDialogComponent(new DialogComponentStringSelection(inputtype,"Type of Input:","Multi-Sample VCFs (Trios)","Single-Sample VCFs"));
    	
    	addDialogComponent(new DialogComponentFileChooser(pedfile, "his_id_ped_denovocaller", 0, ".ped"));
    	
    }
}


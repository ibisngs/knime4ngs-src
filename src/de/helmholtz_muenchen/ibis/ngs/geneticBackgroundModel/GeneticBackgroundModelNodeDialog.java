package de.helmholtz_muenchen.ibis.ngs.geneticBackgroundModel;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "GeneticBackgroundModel" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class GeneticBackgroundModelNodeDialog extends DefaultNodeSettingsPane {
	
	private final SettingsModelString m_resolution = new SettingsModelString(GeneticBackgroundModelNodeModel.CFGKEY_RESOLUTION, GeneticBackgroundModelNodeModel.RESOLUTION[0]);
	private final SettingsModelString m_ac = new SettingsModelString(GeneticBackgroundModelNodeModel.CFGKEY_AC,"AC");
	private final SettingsModelString m_an = new SettingsModelString(GeneticBackgroundModelNodeModel.CFGKEY_AN,"AN");

    /**
     * New pane for configuring the GeneticBackgroundModel node.
     */
    protected GeneticBackgroundModelNodeDialog() {
    	
    	createNewGroup("Resolution");
    	addDialogComponent(new DialogComponentStringSelection(m_resolution, "Compute frequencies using on the level of",GeneticBackgroundModelNodeModel.RESOLUTION));
    	
    	createNewGroup("Field IDs indicating the AC and AN to be used");
    	addDialogComponent(new DialogComponentString(m_ac, "Allele Count (AC)"));
    	addDialogComponent(new DialogComponentString(m_an, "Allele Number (AN)"));
    	
    }
}


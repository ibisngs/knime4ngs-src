package de.helmholtz_muenchen.ibis.ngs.gatkraw;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

/**
 * <code>NodeDialog</code> for the "GATKRaw" Node.
 * GATK node without given command walker
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class GATKRawNodeDialog extends GATKNodeDialog {

    /**
     * New pane for configuring the GATKRaw node.
     */
    protected GATKRawNodeDialog() {

    	
    }

	@Override
	protected void addDialogComponent() {
		final SettingsModelString walker = new SettingsModelString(GATKRawNodeModel.CFGKEY_WALKER,"");
    	final SettingsModelString out = new SettingsModelString(GATKRawNodeModel.CFGKEY_OUT,"");
    	
    	createNewTab("Walker");
    	addDialogComponent(new DialogComponentString(walker, "GATK walker"));
    	
    	createNewGroup("Outfile");
    	addDialogComponent(new DialogComponentFileChooser(out, "hisID_GATK_outfile",0));
		
	}
}

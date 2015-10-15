package de.helmholtz_muenchen.ibis.ngs.depthofcoverage;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

/**
 * <code>NodeDialog</code> for the "DepthOfCoverage" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class DepthOfCoverageNodeDialog extends GATKNodeDialog {

    /**
     * New pane for configuring the DepthOfCoverage node.
     */
	@Override
	protected void addDialogComponent() {

	
		final SettingsModelOptionalString m_extrafilters = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_EXTRAFILTERS,"",false);
		final SettingsModelString m_filesuffix = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_FILESUFFIX,"");
		
		createNewTab("DoC options");
    	
    	addDialogComponent(new DialogComponentString(m_filesuffix, "Outfile Suffix"));
    	addDialogComponent(new DialogComponentOptionalString(m_extrafilters, "GATK Filter Terms (separated by ',')"));
		
	}
}


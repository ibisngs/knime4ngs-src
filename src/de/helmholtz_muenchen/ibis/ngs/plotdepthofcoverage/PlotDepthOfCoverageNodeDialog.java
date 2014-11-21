package de.helmholtz_muenchen.ibis.ngs.plotdepthofcoverage;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "PlotDepthOfCoverage" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class PlotDepthOfCoverageNodeDialog extends DefaultNodeSettingsPane {

	
	private final SettingsModelString m_infolder = new SettingsModelString(PlotDepthOfCoverageNodeModel.CFGKEY_INFOLDER,"");
	private final SettingsModelString m_filesuffix = new SettingsModelString(PlotDepthOfCoverageNodeModel.CFGKEY_FILESUFFIX,"");
	
	
    /**
     * New pane for configuring the PlotDepthOfCoverage node.
     */
    protected PlotDepthOfCoverageNodeDialog() {

    	
    	createNewGroup("Path to Folder with Plot Data");
    	addDialogComponent(new DialogComponentFileChooser(m_infolder, "his_id_GATK_PlotDoC", 0,true,""));
    	addDialogComponent(new DialogComponentString(m_filesuffix, "File Suffix of files to find with Coverage Data)"));
    	
    	
    }
}


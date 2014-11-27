package de.helmholtz_muenchen.ibis.ngs.vcfmerger;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "VCFMerger" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class VCFMergerNodeDialog extends DefaultNodeSettingsPane {

	  private final SettingsModelString INFOLDER = new SettingsModelString(VCFMergerNodeModel.CFGKEY_INFOLDER, "");
	  private final SettingsModelString REGEX = new SettingsModelString(VCFMergerNodeModel.CFGKEY_REGEX, "");
	  private final SettingsModelString OUTFOLDER = new SettingsModelString(VCFMergerNodeModel.CFGKEY_OUTFOLDER, "");

	
	
    /**
     * New pane for configuring the VCFMerger node.
     */
    protected VCFMergerNodeDialog() {
    		
    	createNewGroup("Folder in which the search is performed");
    	addDialogComponent(new DialogComponentFileChooser(INFOLDER, "Folder in which the search is performed", 0,true,""));
    	createNewGroup("Outfolder");
    	addDialogComponent(new DialogComponentFileChooser(OUTFOLDER, "Outfolder", 0,true,""));
    	createNewGroup("");
    	addDialogComponent(new DialogComponentString(REGEX, "File Suffix of VCF Files to merge"));	
    	
    }
}


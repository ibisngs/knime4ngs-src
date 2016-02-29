package de.helmholtz_muenchen.ibis.ngs.gatkCombineVCFs;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


import de.helmholtz_muenchen.ibis.ngs.vcfmerger.VCFMergerNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

/**
 * <code>NodeDialog</code> for the "CombineVCFs" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Kaarin Ahomaa
 */
public class CombineVCFsNodeDialog extends GATKNodeDialog {
	
	 private SettingsModelString GENOTYPEMERGEOPTION;
	 private SettingsModelString OUTFOLDER;


	@Override
	protected void addDialogComponent() {
		
		GENOTYPEMERGEOPTION	= new SettingsModelString(VCFMergerNodeModel.CFGKEY_GENOTYPEMERGEOPTION, "");
		OUTFOLDER = new SettingsModelString(VCFMergerNodeModel.CFGKEY_OUTFOLDER, "");
		
		createNewGroup("GenotypeMergeType");
        addDialogComponent(new DialogComponentStringSelection(GENOTYPEMERGEOPTION, "Genotype Merge Strategy","UNSORTED","UNIQUIFY","REQUIRE_UNIQUE"));
		
        createNewGroup("Folder for output files");
    	addDialogComponent(new DialogComponentFileChooser(OUTFOLDER, "his_id_VEP_OUT", 0, true));
	}
}


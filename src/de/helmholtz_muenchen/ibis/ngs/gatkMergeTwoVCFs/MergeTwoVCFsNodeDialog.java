package de.helmholtz_muenchen.ibis.ngs.gatkMergeTwoVCFs;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

/**
 * <code>NodeDialog</code> for the "MergeTwoVCFs" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Kaarin Ahomaa
 */
public class MergeTwoVCFsNodeDialog extends GATKNodeDialog {
	
	 private SettingsModelString GENOTYPEMERGEOPTION;
	 private SettingsModelString INPUT1;
	 private SettingsModelString INPUT2;
	 private SettingsModelString PRIORITIZE;
	 private SettingsModelString OUTFOLDER;
	 private SettingsModelString FILTEREDRECORDSMERGETYPE;


    @Override
	protected void addDialogComponent() {
		
    	GENOTYPEMERGEOPTION = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_GENOTYPEMERGEOPTION, "");
		INPUT1 = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_INPUT1_TAG, "");
		INPUT2 = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_INPUT2_TAG, "");
		PRIORITIZE = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_PRIORITIZE, "");
		OUTFOLDER = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_OUTFOLDER, "");
		FILTEREDRECORDSMERGETYPE = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_FILTEREDRECORDSMERGETYPE, "");
		
		createNewTab("CombineVariants");
		createNewGroup("GenotypeMergeType");
        addDialogComponent(new DialogComponentStringSelection(GENOTYPEMERGEOPTION, "Genotype Merge Strategy","UNIQUIFY","PRIORITIZE","UNSORTED","REQUIRE_UNIQUE"));
        addDialogComponent(new DialogComponentString(INPUT1, "Input VCF file 1"));
    	addDialogComponent(new DialogComponentString(INPUT2, "Input VCF file 2"));
    	addDialogComponent(new DialogComponentString(PRIORITIZE, "Prioritize input"));
    	
    	createNewGroup("Folder for output files");
    	addDialogComponent(new DialogComponentFileChooser(OUTFOLDER, "his_id_VEP_OUT", 0, true));
    	
    	createNewGroup("FilteredRecordMergeType");
    	addDialogComponent(new DialogComponentStringSelection(FILTEREDRECORDSMERGETYPE, "FilteredRecordMergeType","KEEP_IF_ANY_UNFILTERED","KEEP_IF_ALL_UNFILTERED","KEEP_UNCONDITIONAL"));
	}
}

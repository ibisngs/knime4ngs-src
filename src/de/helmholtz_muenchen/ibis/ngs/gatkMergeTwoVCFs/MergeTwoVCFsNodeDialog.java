package de.helmholtz_muenchen.ibis.ngs.gatkMergeTwoVCFs;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
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


    @Override
	protected void addDialogComponent() {
		
    	GENOTYPEMERGEOPTION = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_GENOTYPEMERGEOPTION, "");
		INPUT1 = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_INPUT1, "");
		INPUT2 = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_INPUT2, "");
		PRIORITIZE = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_PRIORITIZE, "");
		
		createNewGroup("GenotypeMergeType");
        addDialogComponent(new DialogComponentStringSelection(GENOTYPEMERGEOPTION, "Genotype Merge Strategy","UNSORTED","UNIQUIFY","REQUIRE_UNIQUE"));
        addDialogComponent(new DialogComponentString(INPUT1, "Input VCF file 1"));
    	addDialogComponent(new DialogComponentString(INPUT2, "Input VCF file 2"));
    	addDialogComponent(new DialogComponentString(PRIORITIZE, "Prioritize input"));
	}
}

/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
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
		
		setDefaultTabTitle("CombineVariants");
//createNewTab("CombineVariants");
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

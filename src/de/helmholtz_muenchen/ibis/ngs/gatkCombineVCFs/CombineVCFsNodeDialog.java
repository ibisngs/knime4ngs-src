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


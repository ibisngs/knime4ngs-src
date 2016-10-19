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

package de.helmholtz_muenchen.ibis.ngs.summarizelof;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


/**
 * <code>NodeDialog</code> for the "SummarizeLOF" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class SummarizeLOFNodeDialog extends DefaultNodeSettingsPane {

	
	private final SettingsModelString pedfile = new SettingsModelString(SummarizeLOFNodeModel.CFGKEY_PED_FILE,"-");
	private final SettingsModelString vcfin = new SettingsModelString(SummarizeLOFNodeModel.CFGKEY_VCF_INFILE,"-");
	
	
    /**
     * New pane for configuring the SummarizeLOF node.
     */
    protected SummarizeLOFNodeDialog() {
    	
    	createNewGroup("VCF File");
    	addDialogComponent(new DialogComponentFileChooser(vcfin, "his_id_SummarizeLOF_VCFIN", 0, ".vcf"));
    	createNewGroup("PED File");
    	addDialogComponent(new DialogComponentFileChooser(pedfile, "his_id_SummarizeLOF_PED", 0, ".ped")); 	
    }
}


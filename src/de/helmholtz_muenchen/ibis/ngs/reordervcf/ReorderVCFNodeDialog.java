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

package de.helmholtz_muenchen.ibis.ngs.reordervcf;

import javax.swing.JFileChooser;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "ReorderVCF" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class ReorderVCFNodeDialog extends DefaultNodeSettingsPane {

//	private final SettingsModelString INFILE 				= new SettingsModelString(ReorderVCFNodeModel.CFGKEY_INFILE, "---");
	private final SettingsModelString REFERENCE_VCF 		= new SettingsModelString(ReorderVCFNodeModel.CFGKEY_REFERENCE_VCF, "---");
	private final SettingsModelString VCFTOOLS 				= new SettingsModelString(ReorderVCFNodeModel.CFGKEY_VCFTOOLS, "---");


    /**
     * New pane for configuring the ReorderVCF node.
     */
    protected ReorderVCFNodeDialog() {

//    	createNewGroup("Infile");
//    	DialogComponentFileChooser infile= new DialogComponentFileChooser(INFILE, "infile_variant_filter", JFileChooser.OPEN_DIALOG, false, ".vcf");
////    	gatkf.setBorderTitle("Choose File for filtering");
//    	addDialogComponent(infile);
    	createNewGroup("Reference VCF");
    	DialogComponentFileChooser ref_vcf= new DialogComponentFileChooser(REFERENCE_VCF, "infile_variant_filter", JFileChooser.OPEN_DIALOG, false, ".vcf");
//    	gatkf.setBorderTitle("Choose File for filtering");
    	addDialogComponent(ref_vcf);
    	createNewGroup("VCFTools Perl Folder");
    	DialogComponentFileChooser vcftools= new DialogComponentFileChooser(VCFTOOLS, "infile_variant_filter", JFileChooser.OPEN_DIALOG, true, "");
//    	gatkf.setBorderTitle("Choose File for filtering");
    	addDialogComponent(vcftools);
    	
    	
    }
}


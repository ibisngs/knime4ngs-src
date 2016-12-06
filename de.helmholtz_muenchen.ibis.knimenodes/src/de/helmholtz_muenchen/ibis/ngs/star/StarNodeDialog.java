
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
package de.helmholtz_muenchen.ibis.ngs.star;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "Star" Node.
 * STAR aligns RNA-seq reads to a reference genome using uncompressed suffix arrays. 
 * For details, please see paper: 
 * Dobin et al, Bioinformatics 2012; doi: 10.1093/bioinformatics/bts635
 *
 * @author Michael Kluge
 */
public class StarNodeDialog extends HTExecutorNodeDialog {

	private final static String BINARY_NAME = "STAR";
	
    protected StarNodeDialog() {}
    
    public void addToolDialogComponents() {
        final SettingsModelString run_mode				= new SettingsModelString(StarNodeModel.CFGKEY_RUN_MODE, StarNodeModel.DEFAULT_RUN_MODE);
        final SettingsModelBoolean twopass_mode			= new SettingsModelBoolean(StarNodeModel.CFGKEY_TWOPASS, true);
        final SettingsModelString out_folder			= new SettingsModelString(StarNodeModel.CFGKEY_OUTPUT_FOLDER, "");
        final SettingsModelString genome_folder			= new SettingsModelString(StarNodeModel.CFGKEY_GENOME_FOLDER, "");
        final SettingsModelString gtf_file				= new SettingsModelString(StarNodeModel.CFGKEY_GTF_FILE, "");
        final SettingsModelIntegerBounded overhang		= new SettingsModelIntegerBounded(StarNodeModel.CFGKEY_OVERHANG, 100, 1, Integer.MAX_VALUE);
        final SettingsModelIntegerBounded threads 		= new SettingsModelIntegerBounded(StarNodeModel.CFGKEY_THREADS, 4, 1, Integer.MAX_VALUE);
        final SettingsModelOptionalString opt_parameter	= new SettingsModelOptionalString(StarNodeModel.CFGKEY_OPTIONAL_PARA, "",false);
        final SettingsModelString bin_path				= new SettingsModelString(BinaryWrapperNodeModel.CFGKEY_BINARY_PATH, BinaryWrapperNodeModel.DEFAULT_BINARY_PATH);
        
        addPrefPageSetting(bin_path,getNameOfBinary());
        
        DialogComponentFileChooser dcOutputFolder 	= new DialogComponentFileChooser(out_folder, "his_id_OUTPUT_FOLDER", 0, true);
        dcOutputFolder.setBorderTitle("Path to output folder");
        
        DialogComponentFileChooser dcGenomeFolder 	= new DialogComponentFileChooser(genome_folder, "his_id_GENOME_FOLDER", 0, true);
    	dcGenomeFolder.setBorderTitle("Path to genome index");
        
    	DialogComponentFileChooser dcGTFfile 		= new DialogComponentFileChooser(gtf_file, "his_id_gtf_file", 0 , false, ".gtf");
    	dcGTFfile.setBorderTitle("Path to GTF file");
    			
       	DialogComponentStringSelection dcRunMode 	= new DialogComponentStringSelection(run_mode, "runMode:", StarNodeModel.DEFAULT_RUN_MODE, StarNodeModel.ALTERNATIVE_RUN_MODE);

       	// add groups and components
       	createNewGroup("STAR Options");
        addDialogComponent(dcRunMode);
        addDialogComponent(new DialogComponentBoolean(twopass_mode, "Use 2pass mapping?"));
        addDialogComponent(new DialogComponentNumber(threads, "Number of threads", 1));
        addDialogComponent(new DialogComponentOptionalString(opt_parameter,"Optional parameters"));
        
        createNewGroup("Input");
        addDialogComponent(dcGenomeFolder);
        addDialogComponent(dcGTFfile);
        addDialogComponent(new DialogComponentNumber(overhang, "Overhang", 5));

        createNewGroup("Output");
        addDialogComponent(dcOutputFolder);
        
        
        // add change listener to runMode
        run_mode.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				if(StarNodeModel.DEFAULT_RUN_MODE.equals(run_mode.getStringValue())) {
					genome_folder.setEnabled(true);
					twopass_mode.setEnabled(true);
				} else {
					genome_folder.setEnabled(false);
					twopass_mode.setEnabled(false);
				}
			}
        });
    }

	protected String getNameOfBinary() {
		return BINARY_NAME;
	}
}


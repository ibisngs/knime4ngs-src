
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

import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeDialog;

/**
 * <code>NodeDialog</code> for the "Star" Node.
 * STAR aligns RNA-seq reads to a reference genome using uncompressed suffix arrays. 
 * For details, please see paper: 
 * Dobin et al, Bioinformatics 2012; doi: 10.1093/bioinformatics/bts635
 *
 * @author Michael Kluge
 */
public class StarNodeDialog extends BinaryWrapperNodeDialog {

	private final static String BINARY_NAME = "STAR";
	
	
    
    protected StarNodeDialog() {}
    
    public void addToolDialogComponents() {
    	// definition of SettingsModel (all prefixed with SET)
        final SettingsModelString SET_RUN_MODE					= new SettingsModelString(StarNodeModel.CFGKEY_RUN_MODE, StarNodeModel.DEFAULT_RUN_MODE);
        final SettingsModelString SET_OUTPUT_FOLDER				= new SettingsModelString(StarNodeModel.CFGKEY_OUTPUT_FOLDER, StarNodeModel.DEFAULT_OUTPUT_FOLDER);
        final SettingsModelString SET_GENOME_FOLDER				= new SettingsModelString(StarNodeModel.CFGKEY_GENOME_FOLDER, StarNodeModel.DEFAULT_GENOME_FOLDER);
        final SettingsModelOptionalString SET_OPTIONAL_PARA		= new SettingsModelOptionalString(StarNodeModel.CFGKEY_OPTIONAL_PARA, "",false);

        
        // create open file/folder components
        DialogComponentFileChooser dcOutputFolder 	= new DialogComponentFileChooser(SET_OUTPUT_FOLDER, "his_id_OUTPUT_FOLDER", 0, true);
       	DialogComponentFileChooser dcGenomeFolder 	= new DialogComponentFileChooser(SET_GENOME_FOLDER, "his_id_GENOME_FOLDER", 0, true);

       	// create string selection component
       	DialogComponentStringSelection dcRunMode 	= new DialogComponentStringSelection(SET_RUN_MODE, "runMode:", StarNodeModel.DEFAULT_RUN_MODE, StarNodeModel.ALTERNATIVE_RUN_MODE);
       	
       	// set a new title to them
       	dcOutputFolder.setBorderTitle("Path to output folder");
       	dcGenomeFolder.setBorderTitle("Path to genome indexes");
     
       	// add groups and components
       	createNewGroup("STAR Options");
        addDialogComponent(dcRunMode);
        addDialogComponent(new DialogComponentOptionalString(SET_OPTIONAL_PARA,"Optional parameters"));
        
        createNewGroup("Input");
        addDialogComponent(dcGenomeFolder);
        
        createNewGroup("Output");
        addDialogComponent(dcOutputFolder);
        
        // add change listener to runMode
        SET_RUN_MODE.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				if(StarNodeModel.DEFAULT_RUN_MODE.equals(SET_RUN_MODE.getStringValue()))
					SET_GENOME_FOLDER.setEnabled(true);
				else
					SET_GENOME_FOLDER.setEnabled(false);
			}
        });
    }

	@Override
	protected String getNameOfBinary() {
		return BINARY_NAME;
	}
}


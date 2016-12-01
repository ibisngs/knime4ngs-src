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


package de.helmholtz_muenchen.ibis.ngs.matsResultIndexer;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode.BinaryWrapperNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;


/**
 * <code>NodeDialog</code> for the "MatsResultIndexer" Node.
 * 
 * @author Michael Kluge
 */
public class MatsResultIndexerNodeDialog extends HTExecutorNodeDialog {

	private final static String BINARY_NAME = "index_gff";
    
    /**
     * New pane for configuring the MatsResultIndexer node.
     */
    protected MatsResultIndexerNodeDialog() {}
    
    public void addToolDialogComponents() {
    	
    	final SettingsModelString SET_OUTPUT_FILE = new SettingsModelString(MatsResultIndexerNodeModel.CFGKEY_OUTPUT_FILE, MatsResultIndexerNodeModel.DEFAULT_OUTPUT_FOLDER);
        final SettingsModelString SET_INPUT_FILE = new SettingsModelString(MatsResultIndexerNodeModel.CFGKEY_INPUT_FILE, MatsResultIndexerNodeModel.DEFAULT_INPUT_FOLDER);
        final SettingsModelBoolean SET_INCLUDE_NOVEL = new SettingsModelBoolean(MatsResultIndexerNodeModel.CFGKEY_INCLUDE_NOVEL, MatsResultIndexerNodeModel.DEFAULT_INCLUDE_NOVEL);
        final SettingsModelString SET_BINARY_PATH			= new SettingsModelString(BinaryWrapperNodeModel.CFGKEY_BINARY_PATH, BinaryWrapperNodeModel.DEFAULT_BINARY_PATH);

		DialogComponentFileChooser dcInputFile 	= new DialogComponentFileChooser(SET_INPUT_FILE, "his_id_INPUT_FILE_MatsResultIndexer", 0, true);
		DialogComponentFileChooser dcOutputFile 	= new DialogComponentFileChooser(SET_OUTPUT_FILE, "his_id_OUTPUT_FILE_MatsResultIndexer", 0, true);
		DialogComponentBoolean dcIncludeNovel	 	= new DialogComponentBoolean(SET_INCLUDE_NOVEL, "include novel events found by MATS");
		  	
		// set a new title to them
		dcOutputFile.setBorderTitle("Path to output folder");
		dcInputFile.setBorderTitle("Path to input folder");
		
        
		
		createNewTab("MatsResultIndexer");
		    // add groups and components
		createNewGroup("Binary");
		addDialogComponent(new DialogComponentFileChooser(SET_BINARY_PATH, "his_id_BIN", 0, false));
		createNewGroup("Input");
		addDialogComponent(dcInputFile);
		 
		createNewGroup("Output");
		addDialogComponent(dcOutputFile);
		 
		createNewGroup("Further options");
		addDialogComponent(dcIncludeNovel);
    }
    
	protected String getNameOfBinary() {
		return BINARY_NAME;
	}
}


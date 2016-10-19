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


package de.helmholtz_muenchen.ibis.ngs.matsResultPlotter;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "MatsResultPlotter" Node.
 * 
 * @author Michael Kluge
 */
public class MatsResultPlotterNodeDialog extends DefaultNodeSettingsPane {
	
    private final SettingsModelString SET_BINARY_PATH	= new SettingsModelString(MatsResultPlotterNodeModel.CFGKEY_BINARY_PATH, MatsResultPlotterNodeModel.DEFAULT_BINARY_PATH);
    private final SettingsModelString SET_OUTPUT_FILE = new SettingsModelString(MatsResultPlotterNodeModel.CFGKEY_OUTPUT_FILE, MatsResultPlotterNodeModel.DEFAULT_OUTPUT_FOLDER);
    private final SettingsModelString SET_INPUT_FILE = new SettingsModelString(MatsResultPlotterNodeModel.CFGKEY_INPUT_FILE, MatsResultPlotterNodeModel.DEFAULT_INPUT_FOLDER);
    private final SettingsModelString SET_SETTINGS_FILE = new SettingsModelString(MatsResultPlotterNodeModel.CFGKEY_SETTINGS_FILE, MatsResultPlotterNodeModel.DEFAULT_SETTINGS_FILE);
    private final SettingsModelDoubleBounded SET_FDR = new SettingsModelDoubleBounded(MatsResultPlotterNodeModel.CFGKEY_FDR, MatsResultPlotterNodeModel.DEFAULT_FDR, 0, 1);
    private final SettingsModelBoolean SET_INCLUDE_READS = new SettingsModelBoolean(MatsResultPlotterNodeModel.CFGKEY_INCLUDE_READS, MatsResultPlotterNodeModel.DEFAULT_INCLUDE_READS);
    
    /**
     * New pane for configuring the MatsResultPlotter node.
     */
    protected MatsResultPlotterNodeDialog() {
    	DialogComponentFileChooser dcPlotScript		= new DialogComponentFileChooser(SET_BINARY_PATH, "his_id_BINARY_FILE_MatsPlotter", 0, false);
    	DialogComponentFileChooser dcSettingFiles	= new DialogComponentFileChooser(SET_SETTINGS_FILE, "his_id_SETTINGS_FILE_MatsPlotter", 0, false);
    	DialogComponentFileChooser dcInputFile 		= new DialogComponentFileChooser(SET_INPUT_FILE, "his_id_INPUT_FILE_MatsPlotter", 0, true);
		DialogComponentFileChooser dcOutputFile 	= new DialogComponentFileChooser(SET_OUTPUT_FILE, "his_id_OUTPUT_FILE_MatsPlotter", 0, true);
		DialogComponentBoolean dcIncludeReads		= new DialogComponentBoolean(SET_INCLUDE_READS, "include reads and do not only look at junctions");
		DialogComponentNumber dcFDR					= new DialogComponentNumber(SET_FDR, "FDR cutoff", 0.01);
		
		// set a new title to them
		dcOutputFile.setBorderTitle("path to output folder");
		dcInputFile.setBorderTitle("path to input folder (MATS_output)");
		dcPlotScript.setBorderTitle("path to sashimi plot script");
		dcSettingFiles.setBorderTitle("path to settings file");
		  
		    // add groups and components
		createNewGroup("input");
		addDialogComponent(dcPlotScript);
		addDialogComponent(dcSettingFiles);
		addDialogComponent(dcInputFile);
		
		createNewGroup("output");
		addDialogComponent(dcOutputFile);
		 
		createNewGroup("further options");
		addDialogComponent(dcIncludeReads);
		addDialogComponent(dcFDR);
    }
}


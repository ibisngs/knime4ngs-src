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
package de.helmholtz_muenchen.ibis.ngs.featureCountsMerger;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "FeatureCountsMerger" Node.
 * 
 * @author Michael Kluge
 */
public class FeatureCountsMergerNodeDialog extends HTExecutorNodeDialog {

	@Override
	public void addToolDialogComponents() {
		final SettingsModelString SET_OUTPUT_FILE = new SettingsModelString(FeatureCountsMergerNodeModel.CFGKEY_OUTPUT_FILE, FeatureCountsMergerNodeModel.DEFAULT_OUTPUT_FILE);
		final SettingsModelBoolean SET_REMOVE_ENDING = new SettingsModelBoolean(FeatureCountsMergerNodeModel.CFGKEY_REMOVE_ENDING, FeatureCountsMergerNodeModel.DEFAULT_REMOVE_ENDING);
		final SettingsModelBoolean SET_REMOVE_PATH = new SettingsModelBoolean(FeatureCountsMergerNodeModel.CFGKEY_REMOVE_PATH, FeatureCountsMergerNodeModel.DEFAULT_REMOVE_PATH);

		// create open file components
		DialogComponentFileChooser dcOutputFile = new DialogComponentFileChooser(SET_OUTPUT_FILE, "his_id_fcm_OUTPUT_FILE", 0, false);
		DialogComponentBoolean removePath 		= new DialogComponentBoolean(SET_REMOVE_PATH, "Remove the path of the file in the header");
		DialogComponentBoolean removeEnding 	= new DialogComponentBoolean(SET_REMOVE_ENDING, "Remove the ending of the file in the header");
			
		// set a new title to them
		dcOutputFile.setBorderTitle("Path to output file");
			 
		// add groups and components        
		createNewGroup("Output");
		addDialogComponent(dcOutputFile);
			
		createNewGroup("Header format options");
		addDialogComponent(removePath);
		addDialogComponent(removeEnding);
	}
}


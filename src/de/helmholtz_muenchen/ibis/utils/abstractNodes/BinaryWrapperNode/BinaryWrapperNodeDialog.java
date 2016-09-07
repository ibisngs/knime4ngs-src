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
package de.helmholtz_muenchen.ibis.utils.abstractNodes.BinaryWrapperNode;


import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>BinaryWrapperNodeDialog</code> provides some options to set a path to a binary and set additional commands.
 * @author Michael Kluge
 *
 */
public abstract class BinaryWrapperNodeDialog extends HTExecutorNodeDialog {

    // definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_BINARY_PATH			= new SettingsModelString(BinaryWrapperNodeModel.CFGKEY_BINARY_PATH, BinaryWrapperNodeModel.DEFAULT_BINARY_PATH);
//    private final SettingsModelString SET_ADDITIONAL_PARAMETER 	= new SettingsModelString(BinaryWrapperNodeModel.CFGKEY_ADDITIONAL_PARAMETER, BinaryWrapperNodeModel.DEFAULT_ADDITIONAL_PARAMETER);
//    private final SettingsModelString SET_PARAMETER_FILE		= new SettingsModelString(BinaryWrapperNodeModel.CFGKEY_PARAMETER_FILE, BinaryWrapperNodeModel.DEFAULT_PARAMETER_FILE);
   
	protected BinaryWrapperNodeDialog() {
		
		addPrefPageSetting(SET_BINARY_PATH,getNameOfBinary());
		// rename default tab
//		setDefaultTabTitle("General Options");
		
        // create open file/folder components
//        DialogComponentFileChooser dcBinaryPath = new DialogComponentFileChooser(SET_BINARY_PATH, "his_id_BINARY_PATH", 0);
//        DialogComponentFileChooser dcParameterFile 	= new DialogComponentFileChooser(SET_PARAMETER_FILE, "his_id_PARAMETER_FILE", 0);
//        DialogComponentString dcAdditionalParameter = new DialogComponentString(SET_ADDITIONAL_PARAMETER, "Additional parameters");
        
        // set title
//    	dcBinaryPath.setBorderTitle("Path to " + getNameOfBinary() + " binary");
//        dcParameterFile.setBorderTitle("path to parameter file (level 0)");
        
       	// add binary group and components
//        createNewGroup("Binary file path and additional parameters");
//        addDialogComponent(dcBinaryPath);
//        addDialogComponent(dcParameterFile);
//        addDialogComponent(dcAdditionalParameter);
        
        // create new tab for GUI options
//        createNewTab("level 2 options");
	}
	
//	@Override
//	protected void updatePrefs() {
//		if(usePrefPage.getBooleanValue()) {
//	    	String bin_path = IBISKNIMENodesPlugin.getDefault().getToolPathPreference(getNameOfBinary());
//	    	if(bin_path != null && !bin_path.equals("")) {
//	    		SET_BINARY_PATH.setStringValue(bin_path);
//	    		SET_BINARY_PATH.setEnabled(false);
//			} else {
//				SET_BINARY_PATH.setEnabled(true);
//			}
//		} else {
//			SET_BINARY_PATH.setEnabled(true);
//		}
//	}
	
	/**
	 * Should return the name of the binary for the GUI
	 * @return
	 */
	protected abstract String getNameOfBinary();
}

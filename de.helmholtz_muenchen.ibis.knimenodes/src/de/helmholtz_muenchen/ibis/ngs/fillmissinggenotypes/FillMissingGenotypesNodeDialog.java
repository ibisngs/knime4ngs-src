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


package de.helmholtz_muenchen.ibis.ngs.fillmissinggenotypes;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "FillMissingGenotypes" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class FillMissingGenotypesNodeDialog extends DefaultNodeSettingsPane {

	private final SettingsModelString COVERAGEFOLDER 		= new SettingsModelString(FillMissingGenotypesNodeModel.CFGKEY_COVERAGEFILEFOLDER, "");
    private final SettingsModelString FILESUFFIX 			= new SettingsModelString(FillMissingGenotypesNodeModel.CFGKEY_FILESUFFIX, "");

	
	
    /**
     * New pane for configuring the FillMissingGenotypes node.
     */
    protected FillMissingGenotypesNodeDialog() {

    	createNewGroup("Folder in which Coverage files can be found");
    	addDialogComponent(new DialogComponentFileChooser(COVERAGEFOLDER, "Folder in which Coverage files can be found", 0,true,""));
    	createNewGroup("");
    	addDialogComponent(new DialogComponentString(FILESUFFIX, "File Suffix of Coverage Files which are used to fill ./. GTs"));
    	
    }
}


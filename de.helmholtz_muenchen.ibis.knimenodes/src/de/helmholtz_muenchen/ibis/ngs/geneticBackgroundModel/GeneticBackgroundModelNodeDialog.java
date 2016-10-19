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


package de.helmholtz_muenchen.ibis.ngs.geneticBackgroundModel;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "GeneticBackgroundModel" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class GeneticBackgroundModelNodeDialog extends DefaultNodeSettingsPane {
	
	private final SettingsModelString m_resolution = new SettingsModelString(GeneticBackgroundModelNodeModel.CFGKEY_RESOLUTION, GeneticBackgroundModelNodeModel.RESOLUTION[0]);
	private final SettingsModelString m_ac = new SettingsModelString(GeneticBackgroundModelNodeModel.CFGKEY_AC,"AC");
	private final SettingsModelString m_an = new SettingsModelString(GeneticBackgroundModelNodeModel.CFGKEY_AN,"AN");

    /**
     * New pane for configuring the GeneticBackgroundModel node.
     */
    protected GeneticBackgroundModelNodeDialog() {
    	
    	createNewGroup("Resolution");
    	addDialogComponent(new DialogComponentStringSelection(m_resolution, "Compute frequencies using on the level of",GeneticBackgroundModelNodeModel.RESOLUTION));
    	
    	createNewGroup("Field IDs indicating the AC and AN to be used");
    	addDialogComponent(new DialogComponentString(m_ac, "Allele Count (AC)"));
    	addDialogComponent(new DialogComponentString(m_an, "Allele Number (AN)"));
    	
    }
}


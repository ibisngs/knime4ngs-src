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


package de.helmholtz_muenchen.ibis.ngs.vcfSampler;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "VCFSampler" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class VCFSamplerNodeDialog extends DefaultNodeSettingsPane {

	private final SettingsModelInteger ctrls = new SettingsModelInteger(VCFSamplerNodeModel.CFGKEY_CTRLS, 100);
	private final SettingsModelInteger buffer = new SettingsModelInteger(VCFSamplerNodeModel.CFGKEY_BUFFER, 100);
	private final SettingsModelInteger cases = new SettingsModelInteger(VCFSamplerNodeModel.CFGKEY_CASES, 100);
	private final SettingsModelString def = new SettingsModelString(VCFSamplerNodeModel.CFGKEY_DEF,"10/2");
	private final SettingsModelString noise = new SettingsModelString(VCFSamplerNodeModel.CFGKEY_NOISE,"100/1.5");
	
    /**
     * New pane for configuring the VCFSampler node.
     */
    protected VCFSamplerNodeDialog() {

    	addDialogComponent(new DialogComponentNumber(cases, "Number of cases", 100));
    	addDialogComponent(new DialogComponentNumber(ctrls, "Number of controls", 100));
    	addDialogComponent(new DialogComponentString(def, "Define variant frequency changes in cases"));
    	addDialogComponent(new DialogComponentString(noise, "Define sub-population bias"));
    	addDialogComponent(new DialogComponentNumber(buffer, "Buffer size", 50));
    }
}


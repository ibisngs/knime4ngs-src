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
package de.helmholtz_muenchen.ibis.ngs.edgeR;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "EdgeR" Node.
 * 
 * @author Michael Kluge
 */
public class EdgeRNodeDialog extends DefaultNodeSettingsPane {
	
    private final SettingsModelString SET_CORRECTION	= new SettingsModelString(EdgeRNodeModel.CFGKEY_CORRECTION_METHOD, EdgeRNodeModel.DEFAULT_CORRECTION_METHOD);
    private final SettingsModelString SET_NORM_FACTOR	= new SettingsModelString(EdgeRNodeModel.CFGKEY_NORMALIZE_METHOD_FACTOR, EdgeRNodeModel.DEFAULT_NORMALIZE_METHOD_FACTOR);

    /**
	 * Constructor
	 */
    protected EdgeRNodeDialog() {
    	super();
    	// create the components
    	DialogComponentStringSelection cor = new DialogComponentStringSelection(SET_CORRECTION, "P-value correction method", EdgeRNodeModel.CORRECTION_METHODS, false);
    	DialogComponentStringSelection normFac = new DialogComponentStringSelection(SET_NORM_FACTOR, "Method for calculation of normalization factors", EdgeRNodeModel.NORM_FACTORS, false);
    	
    	this.addDialogComponent(normFac);
    	this.addDialogComponent(cor);
    }
}
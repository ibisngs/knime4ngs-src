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
package de.helmholtz_muenchen.ibis.ngs.DESeq;

import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "DESeq" Node.
 * 
 * @author Michael Kluge
 */
public class DESeqNodeDialog extends HTExecutorNodeDialog {
	
	@Override
	public void addToolDialogComponents() {

		final SettingsModelString SET_METHOD	= new SettingsModelString(DESeqNodeModel.CFGKEY_METHOD, DESeqNodeModel.DEFAULT_METHOD);
		final SettingsModelString SET_SHEARING	= new SettingsModelString(DESeqNodeModel.CFGKEY_SHEARING, DESeqNodeModel.DEFAULT_SHEARING);
	
    	DialogComponentStringSelection method = new DialogComponentStringSelection(SET_METHOD, "Empirical dispersion calculation", DESeqNodeModel.METHODS, false);
    	DialogComponentStringSelection shearing = new DialogComponentStringSelection(SET_SHEARING, "Sharing mode", DESeqNodeModel.SHEARING, false);
		
    	this.addDialogComponent(method);
    	this.addDialogComponent(shearing);
    }
}


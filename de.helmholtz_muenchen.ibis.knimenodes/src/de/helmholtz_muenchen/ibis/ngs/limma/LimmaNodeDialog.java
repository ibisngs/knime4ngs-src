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
package de.helmholtz_muenchen.ibis.ngs.limma;

import org.knime.core.data.def.StringCell;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "Limma" Node.
 * 
 * @author Michael Kluge
 */
public class LimmaNodeDialog extends HTExecutorNodeDialog {
	
	public LimmaNodeDialog() {
		super(StringCell.TYPE.getPreferredValueClass(), 1);
	}
	
	@Override
	public void addToolDialogComponents() {
	    final SettingsModelString SET_CORRECTION	= new SettingsModelString(LimmaNodeModel.CFGKEY_CORRECTION_METHOD, LimmaNodeModel.DEFAULT_CORRECTION_METHOD);
	    final SettingsModelString SET_NORM_FACTOR	= new SettingsModelString(LimmaNodeModel.CFGKEY_NORMALIZE_METHOD_FACTOR, LimmaNodeModel.DEFAULT_NORMALIZE_METHOD_FACTOR);
	    final SettingsModelString SET_METHOD_CPM	= new SettingsModelString(LimmaNodeModel.CFGKEY_NORMALIZE_METHOD_CPM, LimmaNodeModel.DEFAULT_NORMALIZE_METHOD_CPM);

	    DialogComponentStringSelection cor = new DialogComponentStringSelection(SET_CORRECTION, "P-value correction method", LimmaNodeModel.CORRECTION_METHODS, false);
	    DialogComponentStringSelection normFac = new DialogComponentStringSelection(SET_NORM_FACTOR, "Method for calculation of normalization factors", LimmaNodeModel.NORM_FACTORS, false);
	    DialogComponentStringSelection normCPM = new DialogComponentStringSelection(SET_METHOD_CPM, "Method for CPM normalization", LimmaNodeModel.NORM_CPM, false);
	    	
	    this.addDialogComponent(normFac);
	    this.addDialogComponent(normCPM);
	    this.addDialogComponent(cor);
		
	}
}


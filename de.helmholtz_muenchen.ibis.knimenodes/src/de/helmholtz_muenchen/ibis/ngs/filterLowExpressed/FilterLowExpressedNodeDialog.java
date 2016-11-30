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
package de.helmholtz_muenchen.ibis.ngs.filterLowExpressed;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentButtonGroup;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;


/**
 * <code>NodeDialog</code> for the "FilterLowExpressed" Node.
 * 
 * @author Michael Kluge
 */
public class FilterLowExpressedNodeDialog extends HTExecutorNodeDialog {
	
	@Override
	public void addToolDialogComponents() {
		
		final SettingsModelInteger SET_KEEP_READS	= new SettingsModelInteger(FilterLowExpressedNodeModel.CFGKEY_KEEP_READS, FilterLowExpressedNodeModel.DEFAULT_KEEP_READS);
		final SettingsModelDouble SET_KEEP_FRACTION	= new SettingsModelDouble(FilterLowExpressedNodeModel.CFGKEY_KEEP_FRACTION, FilterLowExpressedNodeModel.DEFAULT_KEEP_FRATION);
	    final SettingsModelBoolean SET_BOTH_SEP		= new SettingsModelBoolean(FilterLowExpressedNodeModel.CFGKEY_BOTH_SEPERATE, FilterLowExpressedNodeModel.DEFAULT_BOTH_SEP);
	    final SettingsModelString SET_MODE 			= new SettingsModelString(FilterLowExpressedNodeModel.CFGKEY_MODE, FilterLowExpressedNodeModel.DEFAULT_MODE);
		
		// create the components
    	DialogComponentNumber keepReads = new DialogComponentNumber(SET_KEEP_READS, "Minimum read number", 1);
    	DialogComponentNumber keepFraction = new DialogComponentNumber(SET_KEEP_FRACTION, "Fraction of samples", 0.01);
    	DialogComponentBoolean bothSep = new DialogComponentBoolean(SET_BOTH_SEP, "Filter both conditions separately");

		DialogComponentButtonGroup mode = new DialogComponentButtonGroup(SET_MODE, false, "Filtering Mode", FilterLowExpressedNodeModel.DEFAULT_MODE, "Per sample");

    	this.addDialogComponent(mode);
    	this.addDialogComponent(keepReads);
    	this.addDialogComponent(keepFraction);
    	this.addDialogComponent(bothSep);
    	
		// check for changes
    	SET_MODE.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				if(SET_MODE.getStringValue().equals(FilterLowExpressedNodeModel.DEFAULT_MODE)) {
					SET_KEEP_FRACTION.setEnabled(false);
				}
				else {
					SET_KEEP_FRACTION.setEnabled(true);
				}
			}
        });
		
	}
}


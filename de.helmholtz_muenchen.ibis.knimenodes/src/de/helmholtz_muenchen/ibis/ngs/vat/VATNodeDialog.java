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

package de.helmholtz_muenchen.ibis.ngs.vat;

import javax.swing.JFileChooser;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "VAT" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Marie-Sophie Friedl
 */
public class VATNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring VAT node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
	
	final SettingsModelString VATFolder = new SettingsModelString(VATNodeModel.CFGKEY_VAT_FOLDER, VATNodeModel.DEF_VAT_FOLDER);
	final SettingsModelString intervals = new SettingsModelString(VATNodeModel.CFGKEY_INTERVALS, VATNodeModel.DEF_INTERVALS);
	final SettingsModelString transcripts = new SettingsModelString(VATNodeModel.CFGKEY_TRANSCRIPTS, VATNodeModel.DEF_TRANSCRIPTS);

	
    protected VATNodeDialog() {
        super();
        
        createNewGroup("Folder containing VAT executables");
        addDialogComponent(new DialogComponentFileChooser(VATFolder, "vatfolder", JFileChooser.OPEN_DIALOG, true));
        createNewGroup("Path to gene intervals");
        addDialogComponent(new DialogComponentFileChooser(intervals, "geneint", JFileChooser.OPEN_DIALOG, false, ".intervals|.interval"));
        createNewGroup("Path to transcript sequences");
        addDialogComponent(new DialogComponentFileChooser(transcripts, "transcripts", JFileChooser.OPEN_DIALOG, false, ".fa|.fasta"));

                    
    }
}


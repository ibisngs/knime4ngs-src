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
package de.helmholtz_muenchen.ibis.ngs.fastqc_v2;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.ngs.Bcl2FastQ.Bcl2FastQNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "FastQC_v2" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Paul Hager
 */
public class FastQC_v2NodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the FastQC_v2 node.
     */
    protected FastQC_v2NodeDialog() {

    }

	@Override
	public void addToolDialogComponents() {
		//options tab
		final SettingsModelString outfolder = new SettingsModelString(FastQC_v2NodeModel.CFGKEY_OUTFOLDER,"");
		final SettingsModelString fastqc = new SettingsModelString(FastQC_v2NodeModel.CFGKEY_FASTQC,"");
		final SettingsModelIntegerBounded threads = new SettingsModelIntegerBounded(FastQC_v2NodeModel.CFGKEY_THREADS,4, 1, Integer.MAX_VALUE);
		final SettingsModelString additionalOptions = new SettingsModelString(FastQC_v2NodeModel.CFGKEY_ADDITIONAL_OPTIONS,"");
		
		createNewGroup("Folder for output files");
		addDialogComponent(new DialogComponentFileChooser(outfolder, "his_id_FastQC_v2_OUT", 0, true));
		
		createNewGroup("Number of threads");
		addDialogComponent(new DialogComponentNumber(threads, "Number of threads", 1));
		
		createNewGroup("Additional options");
		addDialogComponent(new DialogComponentString(additionalOptions, "Additional FastQC options"));
		
		addPrefPageSetting(fastqc, IBISKNIMENodesPlugin.FASTQC);
	}
}


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
package de.helmholtz_muenchen.ibis.ngs.vepfilter;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentButton;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringListSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;


/**
 * <code>NodeDialog</code> for the "VEPFilter" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class VEPFilterNodeDialog extends HTExecutorNodeDialog {
	
	private DialogComponentStringListSelection DC_TERM_DISPLAY;
	private SettingsModelBoolean my_overwrite;
    
    public void addToolDialogComponents() {
    	
    	final SettingsModelString vep_script = new SettingsModelString(VEPFilterNodeModel.CFGKEY_VEP_SCRIPT,"");
    	final SettingsModelOptionalString filter = new SettingsModelOptionalString(VEPFilterNodeModel.CFGKEY_FILTER,"",false);
    	final SettingsModelString outfolder = new SettingsModelString(VEPFilterNodeModel.CFGKEY_OUTFOLDER,"");
    	final SettingsModelBoolean overwrite = new SettingsModelBoolean(VEPFilterNodeModel.CFGKEY_OVERWRITE,false);
    	my_overwrite = overwrite;
    	
    	final SettingsModelStringArray chosen_terms = new SettingsModelStringArray(VEPFilterNodeModel.CFGKEY_TERM_LIST, new String[]{VEPFilterNodeModel.DEFAULT_TERM});
    	
    	DC_TERM_DISPLAY = new DialogComponentStringListSelection(chosen_terms, "", VEPFilterNodeModel.SO_TERMS);

    	DialogComponentButton CHOOSE_LOF_TERMS = new DialogComponentButton("Choose LOF terms");
    	DialogComponentOptionalString DC_FILTER = new DialogComponentOptionalString(filter,"Conditions");
    	
    	addPrefPageSetting(vep_script, IBISKNIMENodesPlugin.VEP_FILTER);
    	
    	createNewGroup("Filter by Consequence");
		DC_TERM_DISPLAY.setVisibleRowCount(20);
		addDialogComponent(DC_TERM_DISPLAY);
		
		CHOOSE_LOF_TERMS.setToolTipText("Default SO terms comprise variants that are considered as LoF variants.");
		addDialogComponent(CHOOSE_LOF_TERMS);
		
		createNewGroup("Further filtering conditions");
		DC_FILTER.setToolTipText("Define additional comma separated filters,e.g.\"LoF is HC,ExAC_AF>0.05");
		addDialogComponent(DC_FILTER);

		
		CHOOSE_LOF_TERMS.addActionListener(new ActionListener () {
			@Override
			public void actionPerformed(ActionEvent e) {
				chosen_terms.setStringArrayValue(VEPFilterNodeModel.LOF_TERMS);
			}
		});
	
		createNewGroup("Folder for output files");
    	addDialogComponent(new DialogComponentFileChooser(outfolder, "his_id_VEP_Filtered_OUT", 0, true));
    	addDialogComponent(new DialogComponentBoolean(overwrite, "Overwrite, if output files exist?"));
    }
    
    public void onOpen() {
    	super.onOpen();
    	boolean use_hte = IBISKNIMENodesPlugin.getBooleanPreference(IBISKNIMENodesPlugin.USE_HTE);
    	
    	if(use_hte) {
    		my_overwrite.setBooleanValue(true);
    		my_overwrite.setEnabled(false);
    	} else {
//    		my_overwrite.setBooleanValue(false);
    		my_overwrite.setEnabled(true);
    	}
    }
}
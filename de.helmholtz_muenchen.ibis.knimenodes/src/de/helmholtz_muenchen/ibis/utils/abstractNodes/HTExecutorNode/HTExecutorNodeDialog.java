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
package de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode;

import java.awt.Component;
import java.awt.Dimension;
import java.util.LinkedHashMap;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModel;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;

/**
 * <code>NodeDialog</code> for the "HTExecutorNode" Node.
 * 
 * 
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more
 * complex dialog please derive directly from
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public abstract class HTExecutorNodeDialog extends DefaultNodeSettingsPane {

	/**
	 * New pane for configuring the HTExecutorNode node.
	*/
	
	private final SettingsModelBoolean usePrefPage = new SettingsModelBoolean(HTExecutorNodeModel.CFGKEY_USE_PREF,true);
	private final SettingsModelIntegerBounded threshold = new SettingsModelIntegerBounded(HTExecutorNodeModel.CFGKEY_DEFAULT_THRESHOLD, HTExecutorNodeModel.DEFAULT_THRESHOLD, 1, Integer.MAX_VALUE);
	
	protected final LinkedHashMap<SettingsModelString, String> model2pref = new LinkedHashMap<>();
//	private boolean firstOpened;
	
	protected HTExecutorNodeDialog() {
		
		addToolDialogComponents();
		
		if(optionsTabIsEmpty()) {
			setDefaultTabTitle("Preference page");
		} else {
			createNewTab("Preference page");
		}
		
		addDialogComponent(new DialogComponentBoolean(usePrefPage,"Use values from KNIME4NGS preference page?"));
		
		DialogComponentNumber dcn = new DialogComponentNumber(threshold, "HTE threshold", HTExecutorNodeModel.DEFAULT_THRESHOLD);
		addDialogComponent(dcn);
		threshold.setEnabled(false);
		
		usePrefPage.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				threshold.setEnabled(!usePrefPage.getBooleanValue());
				if(usePrefPage.getBooleanValue()) {
					String n = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.THRESHOLD);
			    	threshold.setIntValue(Integer.parseInt(n));
				}
				updatePrefs();
			}
    	});
		
		DialogComponentFileChooser dcfc;
		for(SettingsModelString sm: model2pref.keySet()) {
			createNewGroup("Path to "+model2pref.get(sm));
			dcfc = new DialogComponentFileChooser(sm,model2pref.get(sm),0);
			addDialogComponent(dcfc);
		}
		
//		firstOpened = false;
		
	}
	
    public void onOpen() {
//    	if(!firstOpened) {
//    		DialogComponentFileChooser dcfc;
//    		for(SettingsModelString sm: model2pref.keySet()) {
//    			createNewGroup("Path to "+model2pref.get(sm));
//    			dcfc = new DialogComponentFileChooser(sm,model2pref.get(sm),0);
//    			addDialogComponent(dcfc);
//    		}
//    		firstOpened = true;
//    	}
    	
    	if(usePrefPage.getBooleanValue()) {
			String n = IBISKNIMENodesPlugin.getStringPreference(IBISKNIMENodesPlugin.THRESHOLD);
	    	threshold.setIntValue(Integer.parseInt(n));
	    	threshold.setEnabled(false);
		}
    	updatePrefs();
    }
    
    private void updatePrefs() {
    	for(SettingsModel sm: model2pref.keySet()) {
			sm.setEnabled(true);
		}
    	
    	String prefValue;
    	if(usePrefPage.getBooleanValue()) {
    		for(SettingsModelString sm: model2pref.keySet()) {
    			prefValue = IBISKNIMENodesPlugin.getStringPreference(model2pref.get(sm));
    			if(prefValue != null && !prefValue.equals("")) {
    	    		sm.setStringValue(prefValue);
    	    		sm.setEnabled(false);
    	    	}
    		}
		}
    	
    	for(SettingsModelString sm: model2pref.keySet()) {
    		if(sm.getStringValue().equals("")) {
    			sm.setStringValue("Path is required!");
    		}
    	}
    }
    
    public void addPrefPageSetting(SettingsModelString sms, String v) {
    	this.model2pref.put(sms, v);
    }
    
    public abstract void addToolDialogComponents();
    
    public boolean optionsTabIsEmpty() {
    	Component c = this.getTab("Options");
    	if(c == null) {
    		return false;
    	} 
    	Dimension d = c.getPreferredSize();
    	return (d.getWidth()== 0 && d.getHeight() == 0);
    } 
}

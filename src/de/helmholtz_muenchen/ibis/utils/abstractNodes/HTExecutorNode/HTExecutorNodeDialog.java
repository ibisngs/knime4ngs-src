package de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode;

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
	private boolean firstOpened;
	
	protected HTExecutorNodeDialog() {
		
		addToolDialogComponents();
		
		createNewTab("Preference page");
		
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
		
		firstOpened = false;
		
	}
	
    public void onOpen() {
    	if(!firstOpened) {
    		DialogComponentFileChooser dcfc;
    		for(SettingsModelString sm: model2pref.keySet()) {
    			createNewGroup("Path to "+model2pref.get(sm));
    			dcfc = new DialogComponentFileChooser(sm,model2pref.get(sm),0);
    			addDialogComponent(dcfc);
    		}
    		firstOpened = true;
    	}
    	
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
    }
    
    public void addPrefPageSetting(SettingsModelString sms, String v) {
    	this.model2pref.put(sms, v);
    }
    
    public abstract void addToolDialogComponents();
    
}

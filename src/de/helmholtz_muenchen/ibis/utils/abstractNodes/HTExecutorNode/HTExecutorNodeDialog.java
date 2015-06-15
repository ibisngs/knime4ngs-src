package de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;

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
	private final SettingsModelInteger threshold = new SettingsModelInteger(HTExecutorNodeModel.CFGKEY_DEFAULT_THRESHOLD, HTExecutorNodeModel.DEFAULT_THRESHOLD);
	
	protected HTExecutorNodeDialog() {
		super();

		createNewGroup("HTE");
		
		addDialogComponent(new DialogComponentBoolean(usePrefPage,"Use vales from KNIME4NGS preference page?"));
		
		DialogComponentNumber dcn = new DialogComponentNumber(threshold, "Threshold", HTExecutorNodeModel.DEFAULT_THRESHOLD);
		addDialogComponent(dcn);
		threshold.setEnabled(false);
		usePrefPage.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				threshold.setEnabled(!usePrefPage.getBooleanValue());
				if(usePrefPage.getBooleanValue()) {
					String n = IBISKNIMENodesPlugin.getDefault().getThresholdPreference();
			    	threshold.setIntValue(Integer.parseInt(n));
				}
			}
    	});
	}
	
    public void onOpen() {
    	if(usePrefPage.getBooleanValue()) {
			String n = IBISKNIMENodesPlugin.getDefault().getThresholdPreference();
	    	threshold.setIntValue(Integer.parseInt(n));
	    	threshold.setEnabled(false);
		}
    }
}

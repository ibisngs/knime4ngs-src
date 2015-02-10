package de.helmholtz_muenchen.ibis.htetrigger;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;

/**
 * <code>NodeDialog</code> for the "HTETrigger" Node.
 * 
 * 
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more
 * complex dialog please derive directly from
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class HTETriggerNodeDialog extends DefaultNodeSettingsPane {

	/**
	 * New pane for configuring the HTETrigger node.
	 */
	protected HTETriggerNodeDialog() {
		super();
		addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
				HTETriggerNodeModel.CFGKEY_USE_HTE, HTETriggerNodeModel.USE_HTE),
				"Use HTE?"));
		addDialogComponent(new DialogComponentNumber(new SettingsModelInteger(
				HTETriggerNodeModel.CFGKEY_DEFAULT_THRESHOLD, HTETriggerNodeModel.DEFAULT_THRESHOLD), "Threshold",
				HTETriggerNodeModel.DEFAULT_THRESHOLD));
		addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
				HTETriggerNodeModel.CFGKEY_LOCAL_THRESHOLD, HTETriggerNodeModel.USE_LOCAL_THRESHOLDS),
				"Set individual thresholds?"));
	}
}

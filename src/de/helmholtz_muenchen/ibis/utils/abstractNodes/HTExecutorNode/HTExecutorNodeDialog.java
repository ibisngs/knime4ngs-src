package de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;

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
	protected HTExecutorNodeDialog() {
		super();

		addDialogComponent(new DialogComponentNumber(new SettingsModelInteger(
				HTExecutorNodeModel.CFGKEY_DEFAULT_THRESHOLD, HTExecutorNodeModel.DEFAULT_THRESHOLD), "Threshold",
				HTExecutorNodeModel.DEFAULT_THRESHOLD));

	}
}

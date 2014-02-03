package de.helmholtz_muenchen.ibis.CountMissing;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;

/**
 * <code>NodeDialog</code> for the "CountMissing" Node.
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jonas Zierer
 */
public class CountMissingNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring RNodeTemplate node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
	protected CountMissingNodeDialog() {
        super();

        addDialogComponent(new DialogComponentBoolean(
        		new SettingsModelBoolean(CountMissingNodeModel.CFGKEY_PERCENT, false), "In Percent")
        );


    }
}


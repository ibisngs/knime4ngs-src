package de.helmholtz_muenchen.ibis.misc.sleepnode;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;

/**
 * <code>NodeDialog</code> for the "SleepNode" Node.
 * Node sleeps for the given time
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jonas Zierer
 */
public class SleepNodeNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring SleepNode node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
    protected SleepNodeNodeDialog() {
        super();
        
        addDialogComponent(new DialogComponentNumber(
                new SettingsModelIntegerBounded(
                    SleepNodeNodeModel.CFGKEY_TIME,
                    SleepNodeNodeModel.DEFAULT_TIME,
                    1, 10000000),
                    "Time to sleep:", /*step*/ 1, /*componentwidth*/ 5));
                    
    }
}


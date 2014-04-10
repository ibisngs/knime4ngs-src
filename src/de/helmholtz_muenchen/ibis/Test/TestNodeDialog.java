package de.helmholtz_muenchen.ibis.Test;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;

/**
 * <code>NodeDialog</code> for the "Test" Node.
 * Test node for parallel execution / cluster problem.
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Michael Kluge
 */
public class TestNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring Test node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
    protected TestNodeDialog() {
        super();
                    
    }
}


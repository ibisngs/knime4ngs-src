package de.helmholtz_muenchen.ibis.plotting.barplot;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RPlottingNodeDialog;

/**
 * <code>NodeDialog</code> for the "Boxplot" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jonas Zierer
 */
public class BarplotNodeDialog extends RPlottingNodeDialog {


    /**
     * New pane for configuring the Boxplot node.
     */
	protected BarplotNodeDialog() {
        super();        
        
        this.createNewTabAt("Barplot Options", 0);
  
    }
}


package de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants;


import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

/**
 * <code>NodeDialog</code> for the "GATKSelectVariants" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public abstract class SelectVariantsNodeDialog extends GATKNodeDialog {


    
    /**
     * New pane for configuring the GATKSelectVariants node.
     */

    
    protected void addDialogComponent() {
    	
//    	createNewTab("SelectVariants");
        addFilterDialogComponent();
    }
    
    protected SettingsModelString getFILTERSTRINGModel(){
    	return new SettingsModelString(SelectVariantsNodeModel.CFGKEY_FILTERSTRING, "");
    }
    
    /**
     * Implement specific filter gui
     */
    protected abstract void addFilterDialogComponent();
    
}


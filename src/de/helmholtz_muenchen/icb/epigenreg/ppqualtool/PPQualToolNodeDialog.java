package de.helmholtz_muenchen.icb.epigenreg.ppqualtool;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "PPQualTool" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class PPQualToolNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring PPQualTool node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
    protected PPQualToolNodeDialog() {
        super();
        
        addDialogComponent(new DialogComponentFileChooser(
                new SettingsModelString(
                    PPQualToolNodeModel.CFGKEY_OUTFILE,
                    PPQualToolNodeModel.DEFAULT_OUTFILE), "dubios", 0 , ""));
                    
    }
}


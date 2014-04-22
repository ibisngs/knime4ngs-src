package de.helmholtz_muenchen.ibis.misc.adjustPvalues;

import org.knime.core.data.DoubleValue;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnNameSelection;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "RNodeTemplate" Node.
 * This is a template for a node, which executes a R-Script.
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jonas Zierer
 */
public class PValAdjustNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring RNodeTemplate node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
    
	@SuppressWarnings("unchecked")
	protected PValAdjustNodeDialog() {
        super();
        
        
        this.createNewGroup("Correction Parameters");
        // PValue Column
        addDialogComponent(new DialogComponentColumnNameSelection(
        		new SettingsModelString(PValAdjustNodeModel.CFGKEY_PVAL_COL, "pvalue"),
        		"P-values", 0, DoubleValue.class)
        		);
        // Method Selection
        addDialogComponent(new DialogComponentStringSelection(
        		new SettingsModelString(PValAdjustNodeModel.CFGKEY_METHOD, PValAdjustNodeModel.METHODS[0]), 
        		"Method", PValAdjustNodeModel.METHODS)
        );
        // NUMBER OF TESTS
        addDialogComponent(new DialogComponentNumber(
        		new SettingsModelIntegerBounded(PValAdjustNodeModel.CFGKEY_N, 0, 0, Integer.MAX_VALUE),
        		"Number of Tests", /*step*/ 1, /*componentwidth*/ 5)
        );
        this.closeCurrentGroup();
        
        
        this.createNewGroup("Output Options");
        // Append or Replace Column
        addDialogComponent(new DialogComponentBoolean(
        		new SettingsModelBoolean(PValAdjustNodeModel.CFGKEY_REPLACE, false), "Replace Column")
        );
        addDialogComponent(new DialogComponentString(
        		new SettingsModelString(PValAdjustNodeModel.CFGKEY_COLNAME, "pvalue.corrected"),
        		"New Column Name"
        		)
        );
       this.closeCurrentGroup();     

    }
}


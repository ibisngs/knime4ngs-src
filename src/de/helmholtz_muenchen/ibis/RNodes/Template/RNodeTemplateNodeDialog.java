package de.helmholtz_muenchen.ibis.RNodes.Template;

import org.knime.core.data.DoubleValue;
import org.knime.core.data.IntValue;
import org.knime.core.data.StringValue;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnFilter;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
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
public class RNodeTemplateNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring RNodeTemplate node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
    
	@SuppressWarnings("unchecked")
	protected RNodeTemplateNodeDialog() {
        super();
        
        // INT
        addDialogComponent(new DialogComponentNumber(
        		new SettingsModelIntegerBounded(RNodeTemplateNodeModel.CFGKEY_INT, 5, Integer.MIN_VALUE, Integer.MAX_VALUE),
        		"Counter:", /*step*/ 1, /*componentwidth*/ 5)
        );
        
        // STRING DROPDOWN
        addDialogComponent(new DialogComponentStringSelection(
        		new SettingsModelString(RNodeTemplateNodeModel.CFGKEY_STRING, RNodeTemplateNodeModel.POSSIBLE_STRING_VALUES[0]),
        		"String", /*possible values for dropdown */ RNodeTemplateNodeModel.POSSIBLE_STRING_VALUES)
        		);
        
        // DOUBLE
        addDialogComponent(new DialogComponentNumber(
        		new SettingsModelDouble(RNodeTemplateNodeModel.CFGKEY_DOUBLE, 3.0), "Double Value", /* step size */0.1)
        );
        
        // BOOLEAN
        addDialogComponent(new DialogComponentBoolean(
        		new SettingsModelBoolean(RNodeTemplateNodeModel.CFGKEY_BOOL, false), "Bool Value")
        );
        
        // LIST
        this.createNewGroup("Column name selection");
        addDialogComponent(new DialogComponentColumnFilter(
        		new SettingsModelFilterString(RNodeTemplateNodeModel.CFGKEY_LIST),/* # of input file from which to choose column names*/ 0, true, /* allow columns of the following type */ IntValue.class, DoubleValue.class, StringValue.class)
        );
        this.closeCurrentGroup();

    }
}


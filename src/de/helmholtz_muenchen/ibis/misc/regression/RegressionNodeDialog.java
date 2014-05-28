package de.helmholtz_muenchen.ibis.misc.regression;

import org.knime.core.data.DataValue;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnFilter;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnNameSelection;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "OutlierRemoval" Node.
 * Remove Outliers from Data Columns
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jonas Zierer
 */
public class RegressionNodeDialog extends DefaultNodeSettingsPane {

	/** SETTING MODELS */
	private final SettingsModelString m_target                  = new SettingsModelString(RegressionNodeModel.CFGKEY_TARGET, "");
	private final SettingsModelString m_target_family           = new SettingsModelString(RegressionNodeModel.CFGKEY_TARGET_FAMILY, RegressionNodeModel.COUNFOUNDER_FAMILLIES[0]);
	private final SettingsModelFilterString m_independent       = new SettingsModelFilterString(RegressionNodeModel.CFGKEY_INDEPENDENT);
	private final SettingsModelBoolean m_naaction               = new SettingsModelBoolean(RegressionNodeModel.CFGKEY_NAS, false);

	
    /**
     * New pane for configuring OutlierRemoval node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
    @SuppressWarnings("unchecked")
	protected RegressionNodeDialog() {
        super();
        

		// Variables
		this.createNewGroup("Target Variable");
		this.setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentColumnNameSelection(
				m_target,
				"Select Target Variable",
				0, DataValue.class));
		addDialogComponent(new DialogComponentStringSelection(
				m_target_family,
				"Distribution of target variable"  , RegressionNodeModel.COUNFOUNDER_FAMILLIES)
				);
		this.setHorizontalPlacement(false);
		this.closeCurrentGroup();

		
		this.createNewGroup("Select Independent Variables");
		addDialogComponent(new DialogComponentColumnFilter(
				m_independent,
				0, true, DataValue.class));
		this.closeCurrentGroup();

		this.createNewGroup("Options");
		addDialogComponent(new DialogComponentBoolean(
				m_naaction,
				"Fail on NA"));
		this.closeCurrentGroup();
                    
    }
}


package de.helmholtz_muenchen.ibis.misc.regression.batch;

import org.knime.core.data.DataValue;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnFilter;
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
public class RegressionBatchNodeDialog extends DefaultNodeSettingsPane {

	/** SETTING MODELS */
	private final SettingsModelFilterString m_target            = new SettingsModelFilterString(RegressionBatchNodeModel.CFGKEY_TARGET);
	private final SettingsModelString m_target_family           = new SettingsModelString(RegressionBatchNodeModel.CFGKEY_TARGET_FAMILY, RegressionBatchNodeModel.COUNFOUNDER_FAMILLIES[0]);
	private final SettingsModelFilterString m_independent       = new SettingsModelFilterString(RegressionBatchNodeModel.CFGKEY_INDEPENDENT);
	private final SettingsModelBoolean m_naaction               = new SettingsModelBoolean(RegressionBatchNodeModel.CFGKEY_NAS, false);

	
    /**
     * New pane for configuring OutlierRemoval node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
    @SuppressWarnings("unchecked")
	protected RegressionBatchNodeDialog() {
        super();
        

		// Variables
		this.createNewGroup("Target Variable");
		this.setHorizontalPlacement(false);
		addDialogComponent(new DialogComponentColumnFilter(
				m_target,
				0, false, DataValue.class));
		addDialogComponent(new DialogComponentStringSelection(
				m_target_family,
				"Distribution of target variables"  , RegressionBatchNodeModel.COUNFOUNDER_FAMILLIES)
				);
		this.setHorizontalPlacement(false);
		this.closeCurrentGroup();

		
		this.createNewGroup("Select Independent Variables");
		addDialogComponent(new DialogComponentColumnFilter(
				m_independent,
				0, false, DataValue.class));
		this.closeCurrentGroup();

		this.createNewGroup("Options");
		addDialogComponent(new DialogComponentBoolean(
				m_naaction,
				"Fail on NA"));
		this.closeCurrentGroup();
                    
    }
}


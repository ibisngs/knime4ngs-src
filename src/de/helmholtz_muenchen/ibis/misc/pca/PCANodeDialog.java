package de.helmholtz_muenchen.ibis.misc.pca;

import org.knime.core.data.DataValue;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnFilter;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;

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
public class PCANodeDialog extends DefaultNodeSettingsPane {

	/** SETTING MODELS */
	private final SettingsModelFilterString m_columns   = new SettingsModelFilterString(PCANodeModel.CFGKEY_COLUMNS);
	protected final SettingsModelBoolean m_scale        = new SettingsModelBoolean(PCANodeModel.CFGKEY_SCALE, false);
	protected final SettingsModelBoolean m_center       = new SettingsModelBoolean(PCANodeModel.CFGKEY_CENTER, false);
	private final SettingsModelBoolean m_naaction               = new SettingsModelBoolean(PCANodeModel.CFGKEY_NAS, false);

    /**
     * New pane for configuring OutlierRemoval node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
    @SuppressWarnings("unchecked")
	protected PCANodeDialog() {
        super();

		// Variables
		this.createNewGroup("Select Variables");
		addDialogComponent(new DialogComponentColumnFilter(
				m_columns,
				0, true, DataValue.class));
		this.closeCurrentGroup();

		this.createNewGroup("Data Transformation");
        this.addDialogComponent(new DialogComponentBoolean(
        		m_scale,
                "scale data"));
        this.addDialogComponent(new DialogComponentBoolean(
        		m_center,
                "center data"));
		addDialogComponent(new DialogComponentBoolean(
				m_naaction,
				"Fail on NA"));
		this.closeCurrentGroup();

 
                    
    }
}


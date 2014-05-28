package de.helmholtz_muenchen.ibis.misc.outlierRemoval;

import org.knime.core.data.DataValue;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnFilter;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
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
public class OutlierRemovalNodeDialog extends DefaultNodeSettingsPane {

	/** SETTING MODELS */
	private final SettingsModelString m_method          = new SettingsModelString(OutlierRemovalNodeModel.CFGKEY_METHOD, OutlierRemovalNodeModel.METHODS[0]);
	private final SettingsModelFilterString m_columns   = new SettingsModelFilterString(OutlierRemovalNodeModel.CFGKEY_COLUMNS);
	private final SettingsModelInteger m_devmean_sds    = new SettingsModelInteger(OutlierRemovalNodeModel.CFGKEY_DEVMEAN_SDS, 5);
	
	
    /**
     * New pane for configuring OutlierRemoval node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
    @SuppressWarnings("unchecked")
	protected OutlierRemovalNodeDialog() {
        super();
        
		// method for imputation
		addDialogComponent(new DialogComponentStringSelection(
				m_method,
				"Method for Outlier Detection"  , OutlierRemovalNodeModel.METHODS)
				);

		// Variables
		this.createNewGroup("Select Variables");
		addDialogComponent(new DialogComponentColumnFilter(
				m_columns,
				0, true, DataValue.class));
		this.closeCurrentGroup();
		
		// sds
		this.createNewGroup("Parameters");
		addDialogComponent(new DialogComponentNumber(
				m_devmean_sds, 
				"SDs", 1)
				);
		this.closeCurrentGroup();

 
                    
    }
}


package de.helmholtz_muenchen.ibis.metabolomics.imputing;

import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataValue;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NotConfigurableException;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnFilter;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentStringListSelection;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;

/**
 * <code>NodeDialog</code> for the "Imputer" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jonas Zierer
 */
public class ImputerNodeDialog extends DefaultNodeSettingsPane {

	private DataTableSpec IN_DATA_SPEC;

	@Override
	public void loadAdditionalSettingsFrom(final NodeSettingsRO settings, final DataTableSpec[] specs) throws NotConfigurableException {
		IN_DATA_SPEC = specs[0];

		if(IN_DATA_SPEC==null){
			throw(new NotConfigurableException("No Input data available. Can't configure ImputerNodeDialog!"));
		}
		//System.err.println("VAR TYPES IN DIALOG: " + StringUtils.join(ImputerNodeModel.getVariableTypes(IN_DATA_SPEC), ","));
		String[] selected = null;
		dialog_types.replaceListItems(ImputerNodeModel.getVariableTypes(IN_DATA_SPEC), selected);


	}

	private DialogComponentStringListSelection dialog_types;
	private DialogComponentColumnFilter dialog_columns;

	/** SETTING MODELS */
	private final SettingsModelString m_method          = new SettingsModelString(ImputerNodeModel.CFGKEY_IMPUTING_METHOD, ImputerNodeModel.IMPUTATION_METHODS[0]);
	private final SettingsModelFilterString m_columns   = new SettingsModelFilterString(ImputerNodeModel.CFGKEY_COLUMNS);
	private final SettingsModelStringArray m_varType    = new SettingsModelStringArray(ImputerNodeModel.CFGKEY_VARIABLE_TYPE, new String[]{});
	private final SettingsModelBoolean m_chooseType     = new SettingsModelBoolean(ImputerNodeModel.CFGKEY_CHOOSE_TYPE, false);
	// Parameters for single imputation methods
	private final SettingsModelString m_knn_dist        = new SettingsModelString(ImputerNodeModel.CFGKEY_KNN_DIST, ImputerNodeModel.KNN_DIST_MEAS[0]);
	private final SettingsModelInteger m_knn_k          = new SettingsModelInteger(ImputerNodeModel.CFGKEY_KNN_K, 5);
	private final SettingsModelInteger m_svd_rank       = new SettingsModelInteger(ImputerNodeModel.CFGKEY_RANK, 5);
	private final SettingsModelInteger m_svd_niters     = new SettingsModelInteger(ImputerNodeModel.CFGKEY_NUMITERS, 5);

	@SuppressWarnings("unchecked")
	protected ImputerNodeDialog() {
		super();


		// method for imputation
		addDialogComponent(new DialogComponentStringSelection(
				m_method,
				"Method for Imputation"  , ImputerNodeModel.IMPUTATION_METHODS)
				);

		this.createNewTab("Choose Variables to impute");
		// choose columns directly
		this.createNewGroup("Select Columns");
		dialog_columns = new DialogComponentColumnFilter(
				m_columns,
				0, true, DataValue.class);
		addDialogComponent(dialog_columns);		
		this.closeCurrentGroup();

		// CHOOSE BY TYPE?
		addDialogComponent(new DialogComponentBoolean(
				m_chooseType,
				"Choose by Type"));

		// CHOOSE BY TYPE
		dialog_types = new DialogComponentStringListSelection(
				m_varType,
				"Variable type to impute", new ArrayList<String>(), false, 5);

		addDialogComponent(dialog_types);
		this.closeCurrentGroup();


		this.createNewTab("Imputation Parameters");
		this.createNewGroup("KNN-Imputation");
		addDialogComponent(new DialogComponentNumber(
				m_knn_k, 
				"Knn", 1, 5)
				);
		addDialogComponent(new DialogComponentStringSelection(
				m_knn_dist,
				"Distance Measure", ImputerNodeModel.KNN_DIST_MEAS)
				);

		this.closeCurrentGroup();

		this.createNewGroup("SVD-Imputation");
		addDialogComponent(new DialogComponentNumber(
				m_svd_rank, 
				"Rank k", 1)
				);
		addDialogComponent(new DialogComponentNumber(
				m_svd_niters, 
				"Num. Iterations", 1)
				);
		this.closeCurrentGroup();


		// change listener
		m_method.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				m_knn_k.setEnabled(false);
				m_knn_dist.setEnabled(false);
				m_svd_rank.setEnabled(false);
				m_svd_niters.setEnabled(false);
				if(m_method.getStringValue().equals("knn")){
					m_knn_k.setEnabled(true);
					m_knn_dist.setEnabled(true);
				}else if(m_method.getStringValue().equals("SVD")){
					m_svd_rank.setEnabled(true);
					m_svd_niters.setEnabled(true);
				}
			}
		});

		m_chooseType.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				m_columns.setEnabled(!m_chooseType.getBooleanValue());
				m_varType.setEnabled( m_chooseType.getBooleanValue());
			}
		});
		m_columns.setEnabled(!m_chooseType.getBooleanValue());
		m_varType.setEnabled( m_chooseType.getBooleanValue());
		
		m_varType.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(m_chooseType.getBooleanValue()){
					m_columns.setNewValues(
							ImputerNodeModel.getSelectedColumnsByType(IN_DATA_SPEC, Arrays.asList(m_varType.getStringArrayValue())), 
							ImputerNodeModel.getDeselectedColumnsByType(IN_DATA_SPEC, Arrays.asList(m_varType.getStringArrayValue())), 
							false);
				}
			}
		});
		if(m_chooseType.getBooleanValue()){
			m_columns.setNewValues(
					ImputerNodeModel.getSelectedColumnsByType(IN_DATA_SPEC, Arrays.asList(m_varType.getStringArrayValue())), 
					ImputerNodeModel.getDeselectedColumnsByType(IN_DATA_SPEC, Arrays.asList(m_varType.getStringArrayValue())), 
					false);
		}
	}
}


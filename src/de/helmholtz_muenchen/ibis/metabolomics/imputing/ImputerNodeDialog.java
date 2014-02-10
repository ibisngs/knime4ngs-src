package de.helmholtz_muenchen.ibis.metabolomics.imputing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataValue;
import org.knime.core.data.DoubleValue;
import org.knime.core.data.IntValue;
import org.knime.core.data.StringValue;
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

//	private String[] varTypes = {"1", "2", "3"};// TODO!!!!!
	private DataTableSpec IN_DATA_SPEC;
	
	@Override
	public void loadAdditionalSettingsFrom(final NodeSettingsRO settings, final DataTableSpec[] specs) throws NotConfigurableException {
		System.err.println("LOADING SETTINGS");
		if(specs[0]==null){
			throw(new NotConfigurableException("No Input data available. Can't configure ImputerNodeDialog!"));
		}
		System.err.println("VAR TYPES IN DIALOG: " + StringUtils.join(ImputerNodeModel.getVariableTypes(specs[0]), ","));
		String[] selected = null;
		dialog_types.replaceListItems(ImputerNodeModel.getVariableTypes(specs[0]), selected);
		
		IN_DATA_SPEC = specs[0];
	}

	private DialogComponentStringListSelection dialog_types;
			
	@SuppressWarnings("unchecked")
	protected ImputerNodeDialog() {
		super();

		/** SETTING MODELS */
		final SettingsModelString m_method        = new SettingsModelString(ImputerNodeModel.CFGKEY_IMPUTING_METHOD, ImputerNodeModel.IMPUTATION_METHODS[0]);
		final SettingsModelFilterString m_columns = new SettingsModelFilterString(ImputerNodeModel.CFGKEY_COLUMNS);
		final SettingsModelBoolean m_chooseByType = new SettingsModelBoolean(ImputerNodeModel.CFGKEY_CHOOSE_BY_TYPE, true);
		final SettingsModelStringArray m_varType  = new SettingsModelStringArray(ImputerNodeModel.CFGKEY_VARIABLE_TYPE, new String[]{});
		// Parameters for single imputation methods
		final SettingsModelString m_knn_dist      = new SettingsModelString(ImputerNodeModel.CFGKEY_KNN_DIST, ImputerNodeModel.KNN_DIST_MEAS[0]);
		final SettingsModelInteger m_knn_k        = new SettingsModelInteger(ImputerNodeModel.CFGKEY_KNN_K, 5);
		final SettingsModelInteger m_svd_rank     = new SettingsModelInteger(ImputerNodeModel.CFGKEY_RANK, 5);
		final SettingsModelInteger m_svd_niters   = new SettingsModelInteger(ImputerNodeModel.CFGKEY_NUMITERS, 5);


		// method for imputation
		addDialogComponent(new DialogComponentStringSelection(
				m_method,
				"Method for Imputation"  , ImputerNodeModel.IMPUTATION_METHODS)
				);

		this.createNewTab("Choose Variables to impute");
		// choose columns directly
		this.createNewGroup("Select Columns");
		addDialogComponent(new DialogComponentColumnFilter(
				m_columns,
				0, true, IntValue.class, DoubleValue.class, StringValue.class, DataValue.class)
				);
		this.closeCurrentGroup();
		
		// CHOOSE BY TYPE?
		this.createNewGroup("Choose Columns by Type");
		addDialogComponent(new DialogComponentBoolean(
				m_chooseByType, 
				"Choose columns to impute by Type")
				);
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
		m_chooseByType.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(m_chooseByType.getBooleanValue()){
					m_varType.setEnabled(true);
					m_columns.setEnabled(false);
				}else{
					m_varType.setEnabled(false);
					m_columns.setEnabled(true);
				}
			}
		});
		
		m_varType.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				List<String> columns = ImputerNodeModel.getColumnsByType(IN_DATA_SPEC, Arrays.asList(m_varType.getStringArrayValue()));
				m_columns.setIncludeList(columns);
			}
		});

	}
}


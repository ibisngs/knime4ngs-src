package de.helmholtz_muenchen.ibis.metabolomics.imputing;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * This is the model implementation of Imputer.
 * 
 *
 * @author Jonas Zierer
 */
public class ImputerNodeModel extends RNodeModel {

	static public final String[] IMPUTATION_METHODS = {"random", "min", "max", "mean", "mice", "knn" /*, "SVD", "lm"*/};
	static public final String[] KNN_DIST_MEAS      = {"corr","euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"};
	//public String[] VARIABLE_TYPES  = new String[]{""};


	/** LOGGER */
	private static final NodeLogger LOGGER = NodeLogger.getLogger(ImputerNodeModel.class);

	/** CFG KEYS */
	static public final String CFGKEY_IMPUTING_METHOD       = "Method";
	static public final String CFGKEY_COLUMNS               = "Columns";
	static public final String CFGKEY_CHOOSE_TYPE           = "chose type";
	static public final String CFGKEY_VARIABLE_TYPE         = "Variable Type";
	// Parameters for single imputation methods
	static public final String CFGKEY_KNN_DIST              = "KNN Dist";
	static public final String CFGKEY_KNN_K                 = "KNN K";
	static public final String CFGKEY_RANK                  = "Rank K";
	static public final String CFGKEY_NUMITERS              = "Num of Iterations";


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
	

	/**
	 * Constructor for the node model.
	 */
	protected ImputerNodeModel() {
		super(1, 1,"manipulate" + File.separatorChar + "impute.R", new String[]{"--input"}, new String[]{"--output"});
	}

	/**
	 * {@inheritDoc}
	 * @throws Exception 
	 */
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {

		// check if k < num Rows
		int knn = m_knn_k.getIntValue();
		if(knn > inData[0].getRowCount()){
			LOGGER.info("k is higher than number of rows. Using number of rows -1 (" + (inData[0].getRowCount()-1) + ") instead.");
			knn = inData[0].getRowCount()-1;
		}		

		this.addArgument("--method"   , m_method.getStringValue());
		this.addArgument("--cols"     , m_columns.getIncludeList());
		this.addArgument("--knn"      , knn);
		this.addArgument("--dist"     , m_knn_dist.getStringValue());
		this.addArgument("--rank.k"   , m_svd_rank.getIntValue());
		this.addArgument("--num.iters", m_svd_niters.getIntValue());

		BufferedDataTable[] out = super.execute(inData, exec);
		out[0] = exec.createSpecReplacerTable(out[0], inData[0].getDataTableSpec()); // parse cell types

		return(out);
	}

//	/**
//	 * {@inheritDoc}
//	 */
//	@Override
//	protected void reset() {
//		super.reset();
//	}

	public static List<String> getSelectedColumnsByType(final DataTableSpec inSpecs, final List<String> types){
		ArrayList<String> columns = new ArrayList<String>();
		
		Iterator<DataColumnSpec> it = inSpecs.iterator();
		DataColumnSpec d;
		while(it.hasNext()){
			d = it.next();
			if(types.contains(d.getType().getPreferredValueClass().getName())){
				columns.add(d.getName());
			}
		}
		return(columns);
	}
	public static List<String> getDeselectedColumnsByType(final DataTableSpec inSpecs, final List<String> types){
		ArrayList<String> columns = new ArrayList<String>();
		
		Iterator<DataColumnSpec> it = inSpecs.iterator();
		DataColumnSpec d;
		while(it.hasNext()){
			d = it.next();
			if(!types.contains(d.getType().getPreferredValueClass().getName())){
				columns.add(d.getName());
			}
		}
		return(columns);
	}
	
	public static Set<String> getVariableTypes(final DataTableSpec inSpecs){
		HashMap<String,Boolean> availableTypes = new HashMap<String,Boolean>();

		Iterator<DataColumnSpec> it = inSpecs.iterator();
		DataColumnSpec d;
		while(it.hasNext()){
			d = it.next();
			availableTypes.put(d.getType().getPreferredValueClass().getName(), new Boolean(true));
		}

		return(availableTypes.keySet());
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
		ArrayList<String> newCols       = new ArrayList<String>(Arrays.asList(inSpecs[0].getColumnNames()));
		ArrayList<String> includedCols  = new ArrayList<String>(m_columns.getIncludeList());
		ArrayList<String> excludedCols  = new ArrayList<String>(m_columns.getExcludeList());

		// remove non-existing columns from inlcude-list
		Iterator<String> cit = includedCols.iterator();
		String col;
		while(cit.hasNext()){
			col = cit.next();
			if(! newCols.contains(col)){
				includedCols.remove(col);
			}else{
				newCols.remove(col);
			}
		}
		
		// remove non-existing columns from exclude-list
		cit = excludedCols.iterator();
		while(cit.hasNext()){
			col = cit.next();
			if(! newCols.contains(col)){
				excludedCols.remove(col);
			}else{
				newCols.remove(col);
			}
		}
		
		// add columns which were not there before to exclude list
		cit = newCols.iterator();
		while(cit.hasNext()){
			excludedCols.add(cit.next());
		}
		
		this.m_columns.setNewValues(includedCols, excludedCols, false);
		return(inSpecs);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveSettingsTo(final NodeSettingsWO settings) {
		m_method.saveSettingsTo(settings);
		m_varType.saveSettingsTo(settings);
		m_chooseType.saveSettingsTo(settings);
		m_columns.saveSettingsTo(settings);
		m_knn_dist.saveSettingsTo(settings);
		m_knn_k.saveSettingsTo(settings);
		m_svd_rank.saveSettingsTo(settings);
		m_svd_niters.saveSettingsTo(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
		m_method.loadSettingsFrom(settings);
		m_varType.loadSettingsFrom(settings);
		m_chooseType.loadSettingsFrom(settings);
		m_columns.loadSettingsFrom(settings);
		m_knn_dist.loadSettingsFrom(settings);
		m_knn_k.loadSettingsFrom(settings);
		m_svd_rank.loadSettingsFrom(settings);
		m_svd_niters.loadSettingsFrom(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
		m_method.validateSettings(settings);
		m_varType.validateSettings(settings);
		m_chooseType.validateSettings(settings);
		m_columns.validateSettings(settings);
		m_knn_dist.validateSettings(settings);
		m_knn_k.validateSettings(settings);
		m_svd_rank.validateSettings(settings);
		m_svd_niters.validateSettings(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
		super.loadInternals(internDir, exec);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
		super.saveInternals(internDir, exec);
	}

}


package de.helmholtz_muenchen.ibis.misc.regression.batch;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.IntCell;
import org.knime.core.data.def.StringCell;
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
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * This is the model implementation of OutlierRemover.
 * 
 *
 * @author Jonas Zierer
 */
public class RegressionBatchNodeModel extends RNodeModel {

	public static final String[] COUNFOUNDER_FAMILLIES = new String[]{"binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson"};
    


	/** LOGGER */
	@SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(RegressionBatchNodeModel.class);

	/** CFG KEYS */
	static public final String CFGKEY_TARGET           = "target";
	static public final String CFGKEY_TARGET_FAMILY    = "target family";
	static public final String CFGKEY_INDEPENDENT      = "confounders";
	static public final String CFGKEY_NAS              = "fail on NA";
	

	/** SETTING MODELS */
	private final SettingsModelFilterString m_target            = new SettingsModelFilterString(RegressionBatchNodeModel.CFGKEY_TARGET);
	private final SettingsModelString m_target_family           = new SettingsModelString(RegressionBatchNodeModel.CFGKEY_TARGET_FAMILY, RegressionBatchNodeModel.COUNFOUNDER_FAMILLIES[0]);
	private final SettingsModelFilterString m_independent       = new SettingsModelFilterString(RegressionBatchNodeModel.CFGKEY_INDEPENDENT);
	private final SettingsModelBoolean m_naaction               = new SettingsModelBoolean(RegressionBatchNodeModel.CFGKEY_NAS, false);

	
	/**
	 * Constructor for the node model.
	 */
	protected RegressionBatchNodeModel() {
		super(1, 2,"manipulate" + File.separatorChar + "regressionBatch.R", new String[]{"--input"}, new String[]{"--output", "--residuals"});
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws CanceledExecutionException {
		this.addArgument("--target"     , m_target.getIncludeList());
		this.addArgument("--target.family", m_target_family.getStringValue());
		
		this.addArgument("--confounders", m_independent.getIncludeList());
		this.setFlag("--failOnNA", m_naaction.getBooleanValue());
		
		BufferedDataTable[] out = super.execute(inData, exec);
		out[0] = exec.createSpecReplacerTable(out[0], getResultSpec(inData[0].getDataTableSpec())); // parse cell types
		out[1] = exec.createSpecReplacerTable(out[1], getResidualSpec(inData[0].getDataTableSpec())); // parse cell types
		return(out);
	}


	private DataTableSpec getResultSpec(DataTableSpec inSpec){
		DataColumnSpec[] resultsSpecs = DataTableSpec.createColumnSpecs(new String[]{"target", "p", "n"}, new DataType[]{DataType.getType(StringCell.class), DataType.getType(DoubleCell.class), DataType.getType(IntCell.class)});
		return(new DataTableSpec(resultsSpecs));
	}
	private DataTableSpec getResidualSpec(DataTableSpec inSpec){
		String[] colnames = inSpec.getColumnNames();
		DataType[] types  = new DataType[colnames.length];
		
		for(int i=0; i<colnames.length; i++){
			if(this.m_target.getIncludeList().contains(colnames[i])){
				types[i] = DataType.getType(DoubleCell.class);
			}else{
				types[i] = inSpec.getColumnSpec(colnames[i]).getType();
			}
		}
		DataColumnSpec[] statsSpecs = DataTableSpec.createColumnSpecs(colnames, types);
		return(new DataTableSpec(statsSpecs));
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
		DataTableSpec[] outspecs = new DataTableSpec[]{getResultSpec(inSpecs[0]), getResidualSpec(inSpecs[0])};
		return(outspecs);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveSettingsTo(final NodeSettingsWO settings) {
		m_target.saveSettingsTo(settings);
		m_independent.saveSettingsTo(settings);
		m_target_family.saveSettingsTo(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
		m_target.loadSettingsFrom(settings);
		m_independent.loadSettingsFrom(settings);
		m_target_family.loadSettingsFrom(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
		m_target.validateSettings(settings);
		m_independent.validateSettings(settings);
		m_target_family.validateSettings(settings);
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


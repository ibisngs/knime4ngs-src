package de.helmholtz_muenchen.ibis.misc.regression;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DoubleCell;
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
public class RegressionNodeModel extends RNodeModel {

	public static final String[] COUNFOUNDER_FAMILLIES = new String[]{"binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson"};
    


	/** LOGGER */
	@SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(RegressionNodeModel.class);

	/** CFG KEYS */
	static public final String CFGKEY_TARGET           = "target";
	static public final String CFGKEY_TARGET_FAMILY    = "target family";
	static public final String CFGKEY_INDEPENDENT      = "confounders";
	static public final String CFGKEY_NAS              = "fail on NA";
	

	/** SETTING MODELS */
	private final SettingsModelString m_target                  = new SettingsModelString(RegressionNodeModel.CFGKEY_TARGET, "");
	private final SettingsModelString m_target_family           = new SettingsModelString(RegressionNodeModel.CFGKEY_TARGET_FAMILY, RegressionNodeModel.COUNFOUNDER_FAMILLIES[0]);
	private final SettingsModelFilterString m_independent       = new SettingsModelFilterString(RegressionNodeModel.CFGKEY_INDEPENDENT);
	private final SettingsModelBoolean m_naaction               = new SettingsModelBoolean(RegressionNodeModel.CFGKEY_NAS, false);

	
	/**
	 * Constructor for the node model.
	 */
	protected RegressionNodeModel() {
		super(1, 2,"manipulate" + File.separatorChar + "regression.R", new String[]{"--input"}, new String[]{"--coef", "--residuals"});
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws CanceledExecutionException {
		this.addArgument("--target"     , m_target.getStringValue());
		this.addArgument("--target.family", m_target_family.getStringValue());
		
		this.addArgument("--confounders", m_independent.getIncludeList());
		this.setFlag("--failOnNA", m_naaction.getBooleanValue());
		
		BufferedDataTable[] out = super.execute(inData, exec);
		out[0] = exec.createSpecReplacerTable(out[0], getCoefSpec()); // parse cell types
		out[1] = exec.createSpecReplacerTable(out[1], getResidualSpec()); // parse cell types
		return(out);
	}


	private DataTableSpec getCoefSpec(){
		DataColumnSpec[] statsSpecs = DataTableSpec.createColumnSpecs(new String[]{"Estimate", "Std. Error", "t value", "Pr(>|t|)"}, new DataType[]{DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class)});
		return(new DataTableSpec(statsSpecs));
	}
	private DataTableSpec getResidualSpec(){
		DataColumnSpec[] statsSpecs = DataTableSpec.createColumnSpecs(new String[]{"residuals", "predictors", "fitted", "effects", "weights", "prior.weights"}, new DataType[]{DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class)});
		return(new DataTableSpec(statsSpecs));
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
		DataTableSpec[] outspecs = new DataTableSpec[]{getCoefSpec(), getResidualSpec()};
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


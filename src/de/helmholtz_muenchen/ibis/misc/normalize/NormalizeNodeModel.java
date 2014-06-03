package de.helmholtz_muenchen.ibis.misc.normalize;

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
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * This is the model implementation of OutlierRemover.
 * 
 *
 * @author Jonas Zierer
 */
public class NormalizeNodeModel extends RNodeModel {
	public static final String[] METHODS = new String[]{"quantile normalize", "z-score"};
			
	/** LOGGER */
	@SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(NormalizeNodeModel.class);

	/** CFG KEYS */
	static public final String CFGKEY_METHOD           = "method";
	static public final String CFGKEY_COLUMNS          = "columns";


	/** SETTING MODELS */
	private final SettingsModelString m_method          = new SettingsModelString(NormalizeNodeModel.CFGKEY_METHOD, NormalizeNodeModel.METHODS[0]);
	private final SettingsModelFilterString m_columns   = new SettingsModelFilterString(NormalizeNodeModel.CFGKEY_COLUMNS);


	/**
	 * Constructor for the node model.
	 */
	protected NormalizeNodeModel() {
		super(1, 2,"manipulate" + File.separatorChar + "normalize.R", new String[]{"--input"}, new String[]{"--output", "--stats"});
	}

	/**
	 * {@inheritDoc}
	 * @throws Exception 
	 */
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
		this.addArgument("--cols"     , m_columns.getIncludeList());
		this.addArgument("--method"      , m_method.getStringValue());

		BufferedDataTable[] out = super.execute(inData, exec);
		out[0] = exec.createSpecReplacerTable(out[0], inData[0].getDataTableSpec()); // parse cell types
		out[1] = exec.createSpecReplacerTable(out[1], getStatsSpec()); // parse cell types
		return(out);
	}


    private DataTableSpec getStatsSpec(){
    	DataColumnSpec[] statsSpecs = DataTableSpec.createColumnSpecs(new String[]{"normality.before", "normality.after", "pgain"}, new DataType[]{DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class)});
    	return(new DataTableSpec(statsSpecs));
	}
    
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
		DataTableSpec[] outspecs = new DataTableSpec[]{inSpecs[0], getStatsSpec()};
		return(outspecs);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveSettingsTo(final NodeSettingsWO settings) {
		m_method.saveSettingsTo(settings);
		m_columns.saveSettingsTo(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
		m_method.loadSettingsFrom(settings);
		m_columns.loadSettingsFrom(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
		m_method.validateSettings(settings);
		m_columns.validateSettings(settings);
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


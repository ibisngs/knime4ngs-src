package de.helmholtz_muenchen.ibis.misc.outlierRemoval;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.IntCell;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * This is the model implementation of OutlierRemover.
 * 
 *
 * @author Jonas Zierer
 */
public class OutlierRemovalNodeModel extends RNodeModel {
	public static final String[] METHODS = new String[]{"deviation from mean"};
			
	/** LOGGER */
	@SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(OutlierRemovalNodeModel.class);

	/** CFG KEYS */
	static public final String CFGKEY_METHOD           = "method";
	static public final String CFGKEY_COLUMNS          = "columns";
	static public final String CFGKEY_DEVMEAN_SDS      = "sds";


	/** SETTING MODELS */
	private final SettingsModelString m_method          = new SettingsModelString(OutlierRemovalNodeModel.CFGKEY_METHOD, OutlierRemovalNodeModel.METHODS[0]);
	private final SettingsModelFilterString m_columns   = new SettingsModelFilterString(OutlierRemovalNodeModel.CFGKEY_COLUMNS);
	private final SettingsModelInteger m_devmean_sds    = new SettingsModelInteger(OutlierRemovalNodeModel.CFGKEY_DEVMEAN_SDS, 5);
	

	/**
	 * Constructor for the node model.
	 */
	protected OutlierRemovalNodeModel() {
		super(1, 2,"manipulate" + File.separatorChar + "outlierRemoval.R", new String[]{"--input"}, new String[]{"--output", "--stats"});
	}

	/**
	 * {@inheritDoc}
	 * @throws Exception 
	 */
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
		this.addArgument("--cols"     , m_columns.getIncludeList());
		this.addArgument("--sds"      , m_devmean_sds.getIntValue());

		BufferedDataTable[] out = super.execute(inData, exec);
		out[0] = exec.createSpecReplacerTable(out[0], inData[0].getDataTableSpec()); // parse cell types
		out[1] = exec.createSpecReplacerTable(out[1], getStatsSpec()); // parse cell types
		return(out);
	}


    private DataTableSpec getStatsSpec(){
    	DataColumnSpec[] statsSpecs = DataTableSpec.createColumnSpecs(new String[]{"removed.measurements", "removed.measurements.perc"}, new DataType[]{DataType.getType(IntCell.class), DataType.getType(DoubleCell.class)});
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
		m_devmean_sds.saveSettingsTo(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
		m_method.loadSettingsFrom(settings);
		m_columns.loadSettingsFrom(settings);
		m_devmean_sds.loadSettingsFrom(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
		m_method.validateSettings(settings);
		m_columns.validateSettings(settings);
		m_devmean_sds.validateSettings(settings);
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


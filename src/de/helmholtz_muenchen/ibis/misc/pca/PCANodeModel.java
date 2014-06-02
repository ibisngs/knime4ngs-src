package de.helmholtz_muenchen.ibis.misc.pca;

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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * This is the model implementation of OutlierRemover.
 * 
 *
 * @author Jonas Zierer
 */
public class PCANodeModel extends RNodeModel {
			
	/** LOGGER */
	@SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(PCANodeModel.class);

	/** CFG KEYS */
	static public final String CFGKEY_COLUMNS          = "columns";
	static public final String CFGKEY_SCALE            = "scale";
	static public final String CFGKEY_CENTER           = "center";
	static public final String CFGKEY_NAS              = "fail on NA";

	/** SETTING MODELS */
	private final SettingsModelFilterString m_columns   = new SettingsModelFilterString(PCANodeModel.CFGKEY_COLUMNS);
	protected final SettingsModelBoolean m_scale        = new SettingsModelBoolean(PCANodeModel.CFGKEY_SCALE, false);
	protected final SettingsModelBoolean m_center       = new SettingsModelBoolean(PCANodeModel.CFGKEY_CENTER, false);
	private final SettingsModelBoolean m_naaction               = new SettingsModelBoolean(PCANodeModel.CFGKEY_NAS, false);


	/**
	 * Constructor for the node model.
	 */
	protected PCANodeModel() {
		super(1, 3,"manipulate" + File.separatorChar + "pca.R", new String[]{"--input"}, new String[]{"--output", "--rotation", "--varexplained"});
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
		this.addArgument("--cols"     , m_columns.getIncludeList());
		this.setFlag("--scale", m_scale.getBooleanValue());
		this.setFlag("--center", m_center.getBooleanValue());
		this.setFlag("--failOnNA", m_naaction.getBooleanValue());
		
		BufferedDataTable[] out = super.execute(inData, exec);
		//out[0] = exec.createSpecReplacerTable(out[0], inData[0].getDataTableSpec()); // parse cell types
		out[2] = exec.createSpecReplacerTable(out[2], getVarExplSpec()); // parse cell types
		return(out);
	}


    private DataTableSpec getVarExplSpec(){
    	DataColumnSpec[] statsSpecs = DataTableSpec.createColumnSpecs(new String[]{"num", "lamda", "var.expl", "var.expl.cum"}, new DataType[]{DataType.getType(IntCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class)});
    	return(new DataTableSpec(statsSpecs));
	}
    private DataTableSpec getRotationsSpec(){
    	//DataColumnSpec[] statsSpecs = DataTableSpec.createColumnSpecs(new String[]{"num", "lamda", "var.expl", "var.expl.cum"}, new DataType[]{DataType.getType(IntCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class)});
    	return(null);
	}
    
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
		DataTableSpec[] outspecs = new DataTableSpec[]{null, getRotationsSpec(), getVarExplSpec()};
		return(outspecs);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveSettingsTo(final NodeSettingsWO settings) {
		m_columns.saveSettingsTo(settings);
		m_scale.saveSettingsTo(settings);
		m_center.saveSettingsTo(settings);
		m_naaction.saveSettingsTo(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
		m_columns.loadSettingsFrom(settings);
		m_scale.loadSettingsFrom(settings);
		m_center.loadSettingsFrom(settings);
		m_naaction.loadSettingsFrom(settings);
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
		m_columns.validateSettings(settings);
		m_scale.validateSettings(settings);
		m_center.validateSettings(settings);
		m_naaction.validateSettings(settings);
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


package de.helmholtz_muenchen.ibis.RNodes.Template;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * This is the model implementation of RNodeTemplate.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class RNodeTemplateNodeModel extends RNodeModel {

	public static final String[] POSSIBLE_STRING_VALUES = {"String 1", "String 2", "Another String"};
	
	/** CFG KEYS */
	static final String CFGKEY_INT    = "argInt";
	static final String CFGKEY_STRING = "argString";
	static final String CFGKEY_DOUBLE = "argDouble";
	static final String CFGKEY_LIST   = "argList";
	static final String CFGKEY_BOOL   = "argBoolean";
    
	/** SETTING MODELS */
    private final SettingsModelIntegerBounded m_count  = new SettingsModelIntegerBounded(RNodeTemplateNodeModel.CFGKEY_INT, 5, Integer.MIN_VALUE, Integer.MAX_VALUE);
    private final SettingsModelString         m_string = new SettingsModelString(RNodeTemplateNodeModel.CFGKEY_STRING, RNodeTemplateNodeModel.POSSIBLE_STRING_VALUES[0]);
    private final SettingsModelDouble         m_double = new SettingsModelDouble(RNodeTemplateNodeModel.CFGKEY_DOUBLE, 3.0);
    private final SettingsModelFilterString   m_list   = new SettingsModelFilterString(RNodeTemplateNodeModel.CFGKEY_LIST);
    private final SettingsModelBoolean        m_bool   = new SettingsModelBoolean(RNodeTemplateNodeModel.CFGKEY_BOOL, true);

    
    /**
     * Constructor for the node model.
     */
    protected RNodeTemplateNodeModel() {
        super(1, 2, "utils" + File.separatorChar + "template.R", new String[]{"--input"}, new String[]{"--output", "--output2"});
        // number of input files
        // nuber of output files
        // script path (relative to scripts/R/ directory
        // commandline args for input files
        // commandline args for output files
       
    }

    /**
     * {@inheritDoc}
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws CanceledExecutionException{
    	this.addArgument("--int"   , m_count.getIntValue());
    	this.addArgument("--string", m_string.getStringValue());
    	this.addArgument("--double", m_double.getDoubleValue());
    	this.addArgument("--list"  , m_list.getIncludeList());
    	this.setFlag("--bool", m_bool.getBooleanValue());
    	
		return(super.execute(inData, exec));
	}
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	// TODO
        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
        m_count.saveSettingsTo(settings);
        m_string.saveSettingsTo(settings);
        m_double.saveSettingsTo(settings);
        m_list.saveSettingsTo(settings);
        m_bool.saveSettingsTo(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
        m_count.loadSettingsFrom(settings);
        m_string.loadSettingsFrom(settings);
        m_double.loadSettingsFrom(settings);
        m_list.loadSettingsFrom(settings);
        m_bool.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	m_count.validateSettings(settings);
    	m_string.validateSettings(settings);
        m_double.validateSettings(settings);
        m_list.validateSettings(settings);
        m_bool.validateSettings(settings);

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	super.loadInternals(internDir, exec); // load output from stdout and stderr
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	super.saveInternals(internDir, exec); // save output from stdout and stderr
    }

}


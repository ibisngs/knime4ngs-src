package de.helmholtz_muenchen.ibis.ngs.edgeR;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;

/**
 * This is the model implementation of EdgeR.
 * 
 *
 * @author Michael Kluge
 */
public class EdgeRNodeModel extends RNodeModel {
	
    // keys for SettingsModels
    protected static final String CFGKEY_CORRECTION_METHOD 			= "PvalueCorrection";
    protected static final String CFGKEY_NORMALIZE_METHOD_FACTOR 	= "NormalizeMethodFactor";
    protected static final String CFGKEY_VS							= "VS";
    
    // initial default values for SettingsModels
    protected static final String DEFAULT_CORRECTION_METHOD = "BH";
    protected static final String DEFAULT_NORMALIZE_METHOD_FACTOR = "TMM";
    protected static final String DEFAULT_VS = "";
    
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_CORRECTION	= new SettingsModelString(CFGKEY_CORRECTION_METHOD, DEFAULT_CORRECTION_METHOD);
    private final SettingsModelString SET_NORM_FACTOR	= new SettingsModelString(CFGKEY_NORMALIZE_METHOD_FACTOR, DEFAULT_NORMALIZE_METHOD_FACTOR);
    @SuppressWarnings("unused")
	private final SettingsModelString SET_VS			= new SettingsModelString(CFGKEY_VS, DEFAULT_VS);
    
 
    // the logger instance
	@SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(EdgeRNodeModel.class);
	
	private static final String SCRIPT_PATH = "ngs" + File.separatorChar + "de" + File.separatorChar + "edgeR.R";
	
	// static vars for all availible methods
	protected static final ArrayList<String> CORRECTION_METHODS = new ArrayList<String>();
	protected static final ArrayList<String> NORM_FACTORS = new ArrayList<String>();
	
	static {
		CORRECTION_METHODS.add("BH");
		CORRECTION_METHODS.add("bonferroni");
		CORRECTION_METHODS.add("BY");
		CORRECTION_METHODS.add("fdr");
		CORRECTION_METHODS.add("hochberg");
		CORRECTION_METHODS.add("holm");
		CORRECTION_METHODS.add("hommel");
		CORRECTION_METHODS.add("none");
				
		NORM_FACTORS.add("RLE");
		NORM_FACTORS.add("TMM");
		NORM_FACTORS.add("upperquartile");
		NORM_FACTORS.add("none");
	}

    /**
     * Constructor for the node model.
     */
	protected EdgeRNodeModel() {
		super(2, 1, SCRIPT_PATH, new String[]{"--countTable", "--annotationFile"}, new String[]{"--output"});
	}
	
    /**
     * {@inheritDoc}
     * @throws Exception 
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{    	
		BufferedDataTable[] out = super.execute(inData, exec);
		out[0] = exec.createSpecReplacerTable(out[0], this.getSpec(inData[0].getDataTableSpec())); // parse cell types

		return(out);
	}
    
    /**
     * get specs of table
     * @param inSpec
     * @return
     */
    private DataTableSpec getSpec(DataTableSpec inSpec) {
    	DataColumnSpec[] specs = DataTableSpec.createColumnSpecs(new String[]{"ID", "log2FC", "aveLog2CPM", "PValue", "adj.PValue"}, new DataType[]{DataType.getType(StringCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class)});
    	return(new DataTableSpec(specs));
	}
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	
    	// add values set in GUI
    	this.addArgument("--normFactors", this.SET_NORM_FACTOR.getStringValue());
    	this.addArgument("--correctPvalue", this.SET_CORRECTION.getStringValue());

    	return new DataTableSpec[]{getSpec(inSpecs[0])};
    }

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void reset() {
		super.reset();
	}
	
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
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


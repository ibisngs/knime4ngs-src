package de.helmholtz_muenchen.ibis.ngs.DESeq;

import java.io.File;
import java.util.ArrayList;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;

/**
 * This is the model implementation of DESeq.
 * 
 *
 * @author Michael Kluge
 */
public class DESeqNodeModel extends RNodeModel {
	
    // keys for SettingsModels
    protected static final String CFGKEY_METHOD 	= "method";
    protected static final String CFGKEY_SHEARING 	= "shearing";
    
    // initial default values for SettingsModels
    protected static final String DEFAULT_METHOD = "pooled";
    protected static final String DEFAULT_SHEARING = "maximum";
    
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_METHOD	= new SettingsModelString(CFGKEY_METHOD, DEFAULT_METHOD);
    private final SettingsModelString SET_SHEARING	= new SettingsModelString(CFGKEY_SHEARING, DEFAULT_SHEARING);

    // the logger instance
	@SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(DESeqNodeModel.class);
	
	private static final String SCRIPT_PATH = "ngs" + File.separatorChar + "de" + File.separatorChar + "DESeq.R";
	
	// static vars for all availible methods
	protected static final ArrayList<String> METHODS = new ArrayList<String>();
	protected static final ArrayList<String> SHEARING = new ArrayList<String>();
	
	static {
		METHODS.add("pooled");
		METHODS.add("pooled-CR");
		METHODS.add("per-condition");
		METHODS.add("blind");
		
		SHEARING.add("fit-only");
		SHEARING.add("maximum");
		SHEARING.add("gene-est-only");
	}

    /**
     * Constructor for the node model.
     */
	protected DESeqNodeModel() {
		super(2, 1, SCRIPT_PATH, new String[]{"--countTable", "--annotationFile"}, new String[]{"--output"});
		this.init();
	}
	
	@Override
	public void init() {
		super.init();
		this.addSetting(SET_METHOD);
		this.addSetting(SET_SHEARING);
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
    	DataColumnSpec[] specs = DataTableSpec.createColumnSpecs(new String[]{"ID", "aveLog2CPM", "log2CPM_A", "log2CPM_B", "FC", "log2FC", "PValue", "adj.PValue"}, new DataType[]{DataType.getType(StringCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class)});
    	return(new DataTableSpec(specs));
	}
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {

    	// add values set in GUI
    	this.addArgument("--method", this.SET_METHOD.getStringValue());
    	this.addArgument("--sharingMode", this.SET_SHEARING.getStringValue());
    	return new DataTableSpec[]{getSpec(inSpecs[0])};
    }
}


package de.helmholtz_muenchen.ibis.misc.adjustPvalues;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataTableSpecCreator;
import org.knime.core.data.DataType;
import org.knime.core.data.container.ColumnRearranger;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;


/**
 * This is the model implementation of RNodeTemplate.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class PValAdjustNodeModel extends RNodeModel {

	public static final String[] METHODS = {"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"};
	
	/** CFG KEYS */
	static final String CFGKEY_PVAL_COL = "pvalue column";
	static final String CFGKEY_N    = "number of tests";
	static final String CFGKEY_METHOD = "method";

	static final String CFGKEY_REPLACE   = "replace column";
	static final String CFGKEY_COLNAME   = "new column name";
    
	/** SETTING MODELS */
	private final SettingsModelString         m_pvalcol = new SettingsModelString(PValAdjustNodeModel.CFGKEY_PVAL_COL, "pvalue");
    private final SettingsModelIntegerBounded m_n       = new SettingsModelIntegerBounded(PValAdjustNodeModel.CFGKEY_N, 0, 0, Integer.MAX_VALUE);
    private final SettingsModelString         m_method = new SettingsModelString(PValAdjustNodeModel.CFGKEY_METHOD, "BH");
    
    private final SettingsModelBoolean        m_replace   = new SettingsModelBoolean(PValAdjustNodeModel.CFGKEY_REPLACE, true);
    private final SettingsModelString         m_colname   = new SettingsModelString(PValAdjustNodeModel.CFGKEY_COLNAME, "pvalue.corrected");
    
    /**
     * Constructor for the node model.
     */
    protected PValAdjustNodeModel() {
        super(1, 1, "manipulate" + File.separatorChar + "adjustPvalues.R", new String[]{"--input"}, new String[]{"--output"});
    }

    /**
     * {@inheritDoc}
     * @throws Exception 
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{
    	ColumnRearranger inputRearranger = new ColumnRearranger(inData[0].getDataTableSpec());
    	inputRearranger.keepOnly(m_pvalcol.getStringValue());
    	BufferedDataTable datain = exec.createColumnRearrangeTable(inData[0], inputRearranger, exec);
    	 
    	
    	this.addArgument("--n"     , m_n.getIntValue()==0? inData[0].getRowCount() : m_n.getIntValue()); // if n=0 -> n=number of rows
    	this.addArgument("--method", m_method.getStringValue());
    	this.addArgument("--pvalcol", m_pvalcol.getStringValue());
    	this.addArgument("--outcol" , m_colname.getStringValue());
    	this.addArgument("--format" , "o"); // tell script to return only corrected pvals
    	
    	// Execute script
    	BufferedDataTable[] output = super.execute(new BufferedDataTable[]{datain}, exec);
    	
    	// remove old p-value column (if appropriate)
    	BufferedDataTable dataout;
    	int pcalColIdx = inData[0].getDataTableSpec().findColumnIndex(m_pvalcol.getStringValue()) + 1;
    	if(m_replace.getBooleanValue()){
    		ColumnRearranger removeOldRearranger = new ColumnRearranger(inData[0].getDataTableSpec());
    		removeOldRearranger.remove(m_pvalcol.getStringValue());
    		dataout = exec.createColumnRearrangeTable(inData[0], removeOldRearranger, exec);
    		pcalColIdx = pcalColIdx-1;
    	}else{
    		dataout = inData[0];
    	}
    	
    	// add new column
    	dataout = exec.createJoinedTable(dataout, output[0], exec);
    	
    	
    	// move new column to right place
    	ColumnRearranger moveRearranger = new ColumnRearranger(dataout.getDataTableSpec());
    	moveRearranger.move(m_colname.getStringValue(), pcalColIdx);
    	dataout = exec.createColumnRearrangeTable(dataout, moveRearranger, exec);

    	// cast all types 
    	dataout = exec.createSpecReplacerTable(dataout, this.getOutTableSpec(inData[0].getDataTableSpec())); // parse cell types

		return(new BufferedDataTable[]{dataout});
	}
    
    private DataTableSpec getOutTableSpec(DataTableSpec inSpec){
    	DataTableSpecCreator outSpecCreator = new DataTableSpecCreator(inSpec);
    	outSpecCreator.dropAllColumns();

    	for(String colname: inSpec.getColumnNames()){
    		if(colname.equals(m_pvalcol.getStringValue())){
    			// add p-value column to output if not replaced
    			if(!m_replace.getBooleanValue()){
    				outSpecCreator.addColumns(inSpec.getColumnSpec(colname));
    			}
    			// add new corrected p-value column
    			outSpecCreator.addColumns(DataTableSpec.createColumnSpecs(new String[]{m_colname.getStringValue()}, new DataType[]{DataType.getType(DoubleCell.class)}));
    		}else{
    			// any other column
    			outSpecCreator.addColumns(inSpec.getColumnSpec(colname));
    		}
    	}
    	return(outSpecCreator.createSpec());
	}
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	
        return new DataTableSpec[]{getOutTableSpec(inSpecs[0])};
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
        m_n.saveSettingsTo(settings);
        m_pvalcol.saveSettingsTo(settings);
        m_method.saveSettingsTo(settings);
        m_colname.saveSettingsTo(settings);
        m_replace.saveSettingsTo(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
        m_n.loadSettingsFrom(settings);
        m_pvalcol.loadSettingsFrom(settings);
        m_method.loadSettingsFrom(settings);
        m_colname.loadSettingsFrom(settings);
        m_replace.loadSettingsFrom(settings);
        

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	m_n.validateSettings(settings);
    	m_pvalcol.validateSettings(settings);
        m_method.validateSettings(settings);
        m_colname.validateSettings(settings);
        m_replace.validateSettings(settings);
        
        SettingsModelBoolean m_replace_tmp = m_replace.createCloneWithValidatedValue(settings);
        SettingsModelString  m_colname_tmp = m_colname.createCloneWithValidatedValue(settings);
        SettingsModelString  m_pvalcol_tmp = m_pvalcol.createCloneWithValidatedValue(settings);
        
        if( (!m_replace_tmp.getBooleanValue()) & m_pvalcol_tmp.getStringValue().equals(m_colname_tmp.getStringValue())){
        		throw new InvalidSettingsException("There must not be two columns with the same name.\nEither choose to replace the old p-value column or choose different name for new column!");
        }
        if(m_colname_tmp.getStringValue().equals("")){
        	throw new InvalidSettingsException("Column name of new column must not be empty!");
        }
        

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


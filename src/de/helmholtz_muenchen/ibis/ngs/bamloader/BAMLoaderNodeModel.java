package de.helmholtz_muenchen.ibis.ngs.bamloader;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;


/**
 * This is the model implementation of BAMLoader.
 * 
 *
 * @author 
 */
public class BAMLoaderNodeModel extends NodeModel {

	public static final String CFGKEY_SEQFILE = "seqfile";
	public static final String CFGKEY_BAMFILE = "bamfile";

	private final SettingsModelString m_seqfile = new SettingsModelString(BAMLoaderNodeModel.CFGKEY_SEQFILE,"");
	private final SettingsModelString m_bamfile = new SettingsModelString(BAMLoaderNodeModel.CFGKEY_BAMFILE,"");

    /**
     * Constructor for the node model.
     */
    protected BAMLoaderNodeModel() {
    	
        super(0, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	String path2bamFile = m_bamfile.getStringValue();
    	String path2seqFile = m_seqfile.getStringValue();
    	    	
        DataColumnSpecCreator col2 = new DataColumnSpecCreator("Path2BAMFile", StringCell.TYPE);
        DataColumnSpecCreator col3 = new DataColumnSpecCreator("Path2SEQFile", StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col2.createSpec(),col3.createSpec()};
    	DataTableSpec table = new DataTableSpec(cols);
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	StringCell cl2 = new StringCell(path2bamFile);
    	StringCell cl3 = new StringCell(path2seqFile);
    	DataCell[] c = new DataCell[]{cl2,cl3};
    	DefaultRow r = new DefaultRow("Row0",c);
    	cont.addRowToTable(r);
    	cont.close();
    	BufferedDataTable out = cont.getTable();
    	
    	//Extra FlowVariables
    	pushFlowVariableString("BAMFILE",path2bamFile);
    	pushFlowVariableString("BAMSAMINFILE",path2bamFile);
    	pushFlowVariableString("MPILEUPINFILE", path2bamFile);
    	pushFlowVariableString("Path2seqFile", m_seqfile.getStringValue());
    	
   	
        return new BufferedDataTable[]{out};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	if(m_seqfile.getStringValue().length()>0){
        	if(!FileValidator.checkFastaFormat(m_seqfile.getStringValue())){
                throw new InvalidSettingsException("Reference (genome) sequence file is not in FastA format or does not contain nucleotide sequences!");
        	}
    	}	
    		
    	//Outspecs
    	DataColumnSpecCreator col2 = new DataColumnSpecCreator("Path2BAMFile", StringCell.TYPE);
        DataColumnSpecCreator col3 = new DataColumnSpecCreator("Path2SEQFile", StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col2.createSpec(),col3.createSpec()};
    	DataTableSpec table = new DataTableSpec(cols);
    	
        return new DataTableSpec[]{table};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_seqfile.saveSettingsTo(settings);
    	m_bamfile.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_seqfile.loadSettingsFrom(settings);
    	m_bamfile.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_seqfile.validateSettings(settings);
    	m_bamfile.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    	
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    	
    }

}


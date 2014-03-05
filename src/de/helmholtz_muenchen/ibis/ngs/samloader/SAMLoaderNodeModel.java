package de.helmholtz_muenchen.ibis.ngs.samloader;

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
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;

/**
 * This is the model implementation of SAMLoader.
 * 
 */
public class SAMLoaderNodeModel extends NodeModel {
    
	public static final String CFGKEY_SEQFILE = "seqfile";
	public static final String CFGKEY_SAMFILE = "samfile";

	private final SettingsModelString m_seqfile = new SettingsModelString(SAMLoaderNodeModel.CFGKEY_SEQFILE,"");
	private final SettingsModelString m_samfile = new SettingsModelString(SAMLoaderNodeModel.CFGKEY_SAMFILE,"");
	
    /**
     * Constructor for the node model.
     */
    protected SAMLoaderNodeModel() {
    	
    	super(0, 1);
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	String path2samFile = m_samfile.getStringValue();
    	String path2seqFile = m_seqfile.getStringValue();
    	
    	/**Initialize logfile**/
    	String logfile = path2seqFile.substring(0,path2seqFile.lastIndexOf("/")+1)+"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("SAMLoader"));
    	/**end initializing logfile**/
    	
    	File file1 = new File(path2samFile);
    	File file2 = new File(path2seqFile);
		if(file1.exists() && file2.exists()){
			logBuffer.append("Sequence file and SAM file found.\n");
		}
    	
    	DataColumnSpecCreator col = new DataColumnSpecCreator("Path2SAMFile", StringCell.TYPE);
    	DataColumnSpecCreator col1 = new DataColumnSpecCreator("Sequence file", StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col.createSpec(),col1.createSpec()};
    	DataTableSpec table = new DataTableSpec(cols);
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	StringCell cl1 = new StringCell(path2samFile);
    	StringCell cl2 = new StringCell(path2seqFile);
    	DataCell[] c = new DataCell[]{cl1,cl2};
    	DefaultRow r = new DefaultRow("Row0",c);
    	cont.addRowToTable(r);
    	cont.close();
    	BufferedDataTable out = cont.getTable();
    	
    	pushFlowVariableString("BAMSAMINFILE",path2samFile);
    	
    	logBuffer.append(ShowOutput.getNodeEndTime());
    	ShowOutput.writeLogFile(logBuffer);
    	
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
        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_seqfile.saveSettingsTo(settings);
    	m_samfile.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_seqfile.loadSettingsFrom(settings);
    	m_samfile.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_seqfile.validateSettings(settings);
    	m_samfile.validateSettings(settings);
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

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
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;

/**
 * This is the model implementation of BAMLoader.
 * 
 *
 * @author 
 */
public class BAMLoaderNodeModel extends NodeModel {

	public static final String CFGKEY_SEQFILE = "seqfile";
	public static final String CFGKEY_BAMFILE = "bamfile";
	public static final String CFGKEY_SAMTOOLS = "samtools";

	private final SettingsModelString m_seqfile = new SettingsModelString(BAMLoaderNodeModel.CFGKEY_SEQFILE,"");
	private final SettingsModelString m_bamfile = new SettingsModelString(BAMLoaderNodeModel.CFGKEY_BAMFILE,"");
	private final SettingsModelString m_samtools = new SettingsModelString(BAMLoaderNodeModel.CFGKEY_SAMTOOLS,"");

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
    	
    	String path2samtools = m_samtools.getStringValue();
    	String path2bamFile = m_bamfile.getStringValue();
    	String path2seqFile = m_seqfile.getStringValue();
    	
    	/**Initialize logfile**/
    	String logfile = path2seqFile.substring(0,path2seqFile.lastIndexOf("/")+1)+"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("BAMLoader"));
    	/**end initializing logfile**/
    	
    	File file1 = new File(path2bamFile);
    	File file2 = new File(path2seqFile);
		if(file1.exists() && file2.exists()){
			logBuffer.append("Sequence file and BAM file found.\nVersion of SamTools checked.\n");
		}
    	
    	DataColumnSpecCreator col1 = new DataColumnSpecCreator("Path2SamTools", StringCell.TYPE);
        DataColumnSpecCreator col2 = new DataColumnSpecCreator("Path2BAMFile", StringCell.TYPE);
        DataColumnSpecCreator col3 = new DataColumnSpecCreator("Path2SEQFile", StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec(),col2.createSpec(),col3.createSpec()};
    	DataTableSpec table = new DataTableSpec(cols);
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	StringCell cl1 = new StringCell(path2samtools);
    	StringCell cl2 = new StringCell(path2bamFile);
    	StringCell cl3 = new StringCell(path2seqFile);
    	DataCell[] c = new DataCell[]{cl1,cl2,cl3};
    	DefaultRow r = new DefaultRow("Row0",c);
    	cont.addRowToTable(r);
    	cont.close();
    	BufferedDataTable out = cont.getTable();
    	//Extra FlowVariables
    	pushFlowVariableString("BAMFILE",path2bamFile);
    	pushFlowVariableString("BAMSAMINFILE",path2bamFile);
    	pushFlowVariableString("MPILEUPINFILE", path2bamFile);
    	pushFlowVariableString("Path2SamTools", m_samtools.getStringValue());
    	pushFlowVariableString("Path2seqFile", m_seqfile.getStringValue());
    	
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
    	//Version control
    	if(m_samtools.getStringValue().length()>0){
    		try {
	    		 //Version control
	            if(FileValidator.versionControl(m_samtools.getStringValue(),"SAMTOOLS")==1){
	            	setWarningMessage("WARNING: You are using a newer SAMTOOLS version than "+FileValidator.SAMTOOLS_VERSION +"! This may cause problems!");
	            }else if(FileValidator.versionControl(m_samtools.getStringValue(),"SAMTOOLS")==2){
	            	setWarningMessage("WARNING: You are using an older SAMTOOLS version than "+FileValidator.SAMTOOLS_VERSION +"! This may cause problems!");
	            }else if(FileValidator.versionControl(m_samtools.getStringValue(),"SAMTOOLS")==-1){
	            	setWarningMessage("Your samtools version could not be determined! Correct behaviour can only be ensured for samtools version "+FileValidator.SAMTOOLS_VERSION+".");
	            }
    		} catch (Exception e) {
    			throw new InvalidSettingsException("Specify a valid SAMTOOLS version!");
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
    	m_bamfile.saveSettingsTo(settings);
    	m_samtools.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_seqfile.loadSettingsFrom(settings);
    	m_bamfile.loadSettingsFrom(settings);
    	m_samtools.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_seqfile.validateSettings(settings);
    	m_bamfile.validateSettings(settings);
    	m_samtools.validateSettings(settings);
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


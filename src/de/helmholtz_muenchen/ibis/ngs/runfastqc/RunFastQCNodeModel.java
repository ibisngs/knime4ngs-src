package de.helmholtz_muenchen.ibis.ngs.runfastqc;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
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

import de.helmholtz_muenchen.ibis.utils.BinaryHandler;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;

/**
 * This is the model implementation of RunFastQC.
 * 
 *
 * @author JMB
 */
public class RunFastQCNodeModel extends NodeModel {
	
	public static final String CFGKEY_READSEQFILE = "readseqfile";
	public static final String CFGKEY_READSEQFILE2 = "readseqfile2";
	public static final String CFGKEY_READTYPE = "readType";
	
	private final SettingsModelString m_readseqfile = new SettingsModelString(RunFastQCNodeModel.CFGKEY_READSEQFILE,"");
	private final SettingsModelString m_readseqfile2 = new SettingsModelString(RunFastQCNodeModel.CFGKEY_READSEQFILE2,"");
	private final SettingsModelString m_readType = new SettingsModelString(RunFastQCNodeModel.CFGKEY_READTYPE, "single-end");

	//The Output Col Names
	public static final String OUT_COL1 = "Path2ReadFile1";
	public static final String OUT_COL2 = "Path2ReadFile2";
	
	
    /**
     * Constructor for the node model.
     */
    protected RunFastQCNodeModel() {
    	
        super(0, 1);
        
        m_readseqfile2.setEnabled(false);
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	BinaryHandler.checkToolAvailability("bwa");
    	
    	String readFile1 = m_readseqfile.getStringValue();
    	String readFile2 = "";
    	String readType = m_readType.getStringValue();
    	
    	/**Initialize logfile**/
    	String logfile = readFile1.substring(0,readFile1.lastIndexOf("/")+1)+"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("RunFastQC"));
    	/**end initializing logfile**/
    	
    	/**
    	 * Check Files
    	 */
    	File file1 = new File(readFile1);
    	File file3 = file1;
    	if(m_readseqfile2.isEnabled()) {
    		readFile2 = m_readseqfile2.getStringValue();
    		file3 = new File(readFile2);
    	}
		if(file1.exists() && file3.exists()){
			logBuffer.append("Sequence file(s) found.\n");
		}
		
		/**
		 * Create Output Table
		 */
    	DataColumnSpecCreator col1 = new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE);
        DataColumnSpecCreator col2 = new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE);
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    					new DataColumnSpec[]{
    							col1.createSpec(),col2.createSpec()
    							}));
    	FileCell cl1 = (FileCell) FileCellFactory.create(readFile1);
    	FileCell cl2 = (FileCell) FileCellFactory.create(readFile2);
    	cont.addRowToTable(new DefaultRow("Row0",new FileCell[]{cl1,cl2}));
    	cont.close();
    	BufferedDataTable out = cont.getTable();
    	
    	/**
    	 * Check if Inputformat is bam
    	 */
    	String secFile = "";
    	if(readFile1.substring(readFile1.length()-3,readFile1.length()).equals("bam")) {
    		secFile = "true";
    	} else {
    		secFile = "false";
    	}
    	
    	/**
    	 * Save to FlowVariables
    	 */
    	pushFlowVariableString("isBAM", secFile);
    	pushFlowVariableString("readType", readType);
    	
    	logBuffer.append(ShowOutput.getNodeEndTime());
    	ShowOutput.writeLogFile(logBuffer);
    	
        return new BufferedDataTable[]{out};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
    	//Nothing to do...
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	String seq1 = m_readseqfile.getStringValue();
    	if(seq1.length() > 1) {
	    	if(!seq1.substring(seq1.length()-3,seq1.length()).equals("bam")) {
	    		if(!FileValidator.checkFastaFormat(seq1) && !FileValidator.checkFastqFormat(seq1)){
	                throw new InvalidSettingsException("Read sequences file no. 1 is not in Bam, FastQ or FastA format or does not contain nucleotide sequences!");
	        	}
	    	}
    	}
    	
    	if(m_readseqfile2.isEnabled()&&m_readseqfile2.getStringValue().length()>0) {
    		String seq2 = m_readseqfile2.getStringValue();
    		if(seq2.length() > 1) {
	        	if(!FileValidator.checkFastaFormat(seq2) && !FileValidator.checkFastqFormat(seq2)){
	        		throw new InvalidSettingsException("Read sequences file no. 2 is not in FastQ or FastA format or does not contain nucleotide sequences!");
	        	}
    		}
    	}
    	
        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_readseqfile.saveSettingsTo(settings);
    	m_readseqfile2.saveSettingsTo(settings);
    	m_readType.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_readseqfile.loadSettingsFrom(settings);
    	m_readseqfile2.loadSettingsFrom(settings);
    	m_readType.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_readseqfile.validateSettings(settings);
    	m_readseqfile2.validateSettings(settings);
    	m_readType.validateSettings(settings);
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


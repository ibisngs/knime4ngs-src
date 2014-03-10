package de.helmholtz_muenchen.ibis.ngs.genscan;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.ngs.bowtie2.Bowtie2NodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.QSub;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of Genscan.
 * 
 */
public class GenscanNodeModel extends NodeModel {
    
	public static final String CFGKEY_GENSCANFILE = "genscanfile";
	public static final String CFGKEY_USEMATRIX = "usematrix";
	public static final String CFGKEY_MATRIXFILE = "matrixfile";
	public static final String CFGKEY_VERBOSEOUTPUT = "verboseoutput";
	public static final String CFGKEY_CDS = "cds";
	public static final String CFGKEY_SUBOPT = "subopt";
	public static final String CFGKEY_SUBOPTCUTOFF = "suboptcutoff";
	public static final String CFGKEY_PS = "ps";
	public static final String CFGKEY_PSSCALE = "psscale";

	private final SettingsModelString m_genscanfile = new SettingsModelString(GenscanNodeModel.CFGKEY_GENSCANFILE,"");
	private final SettingsModelBoolean m_usematrix = new SettingsModelBoolean(GenscanNodeModel.CFGKEY_USEMATRIX, true);
	private final SettingsModelString m_matrixfile = new SettingsModelString(GenscanNodeModel.CFGKEY_MATRIXFILE,"");
	private final SettingsModelBoolean m_verboseoutput = new SettingsModelBoolean(GenscanNodeModel.CFGKEY_VERBOSEOUTPUT, true);
	private final SettingsModelBoolean m_cds = new SettingsModelBoolean(GenscanNodeModel.CFGKEY_CDS, true);
	private final SettingsModelBoolean m_subopt = new SettingsModelBoolean(GenscanNodeModel.CFGKEY_SUBOPT, false);
	private final SettingsModelDoubleBounded m_suboptcutoff = new SettingsModelDoubleBounded(GenscanNodeModel.CFGKEY_SUBOPTCUTOFF, 0.1, 0.01, 0.99);
	private final SettingsModelBoolean m_ps = new SettingsModelBoolean(GenscanNodeModel.CFGKEY_PS, true);
	private final SettingsModelIntegerBounded m_psscale = new SettingsModelIntegerBounded(GenscanNodeModel.CFGKEY_PSSCALE, 1, 1, 100);

	private static final NodeLogger LOGGER = NodeLogger.getLogger(GenscanNodeModel.class);
	
	
    /**
     * Constructor for the node model.
     */
    protected GenscanNodeModel() {
        
    	super(1, 0);
        
    	m_matrixfile.setEnabled(false);
    	m_suboptcutoff.setEnabled(false);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	ArrayList<String> command = new ArrayList<String>();
  
    	String path2genscan = m_genscanfile.getStringValue();
    	String model = path2genscan.substring(0,path2genscan.lastIndexOf("/")+1)+"HumanIso.smat";
    	String path2seqfile =  inData[0].iterator().next().getCell(0).toString();
    	String path2outfile = path2seqfile.substring(0,path2seqfile.lastIndexOf("."))+"_genes.txt";
    	String path2outfileps = path2seqfile.substring(0,path2seqfile.lastIndexOf("."))+"_genes.ps";
    	
    	/**Initialize logfile**/
    	String logfile = path2seqfile.substring(0,path2seqfile.lastIndexOf(".")+1)+"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("Genscan"));
    	/**end initializing logfile**/
    	
    	command.add(path2genscan);
    	if(!m_usematrix.getBooleanValue()) {
    		model = m_matrixfile.getStringValue();
    	}
    	command.add(model);
    	command.add(path2seqfile);
    	if(m_verboseoutput.getBooleanValue()) {
    		command.add("-v");
    	}
    	if(m_cds.getBooleanValue()) {
    		command.add("-cds");
    	}
    	if(m_ps.getBooleanValue()) {
    		command.add("-ps " + path2outfileps + " " + m_psscale.getIntValue());
    	}
    	if(m_subopt.getBooleanValue()) {
    		command.add("-subopt " + m_suboptcutoff.getDoubleValue());
    	}
    	
     	/**
     	 * Execute
     	 */
    	// genscan /home/q/quell/tools/GENSCAN/HumanIso.smat sequence.fasta -v -cds -ps out.ps 1 > out.txt
     	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER,path2outfile);
     	logBuffer.append(ShowOutput.getNodeEndTime());
     	ShowOutput.writeLogFile(logBuffer);
        return new BufferedDataTable[]{};
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
    	
        return new DataTableSpec[]{};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_genscanfile.saveSettingsTo(settings);
    	m_usematrix.saveSettingsTo(settings);
    	m_matrixfile.saveSettingsTo(settings);
    	m_verboseoutput.saveSettingsTo(settings);
    	m_cds.saveSettingsTo(settings);
    	m_subopt.saveSettingsTo(settings);
    	m_suboptcutoff.saveSettingsTo(settings);
    	m_ps.saveSettingsTo(settings);
    	m_psscale.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_genscanfile.loadSettingsFrom(settings);
    	m_usematrix.loadSettingsFrom(settings);
    	m_matrixfile.loadSettingsFrom(settings);
    	m_verboseoutput.loadSettingsFrom(settings);
    	m_cds.loadSettingsFrom(settings);
    	m_subopt.loadSettingsFrom(settings);
    	m_suboptcutoff.loadSettingsFrom(settings);
    	m_ps.loadSettingsFrom(settings);
    	m_psscale.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_genscanfile.validateSettings(settings);
    	m_usematrix.validateSettings(settings);
    	m_matrixfile.validateSettings(settings);
    	m_verboseoutput.validateSettings(settings);
    	m_cds.validateSettings(settings);
    	m_subopt.validateSettings(settings);
    	m_suboptcutoff.validateSettings(settings);
    	m_ps.validateSettings(settings);
    	m_psscale.validateSettings(settings);
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


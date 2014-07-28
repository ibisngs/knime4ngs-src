package de.helmholtz_muenchen.ibis.ngs.depthofcoverage;

import java.io.File;
import java.io.IOException;

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
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.ngs.gatkrealignment.GATKRealignmentNodeModel;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;



/**
 * This is the model implementation of DepthOfCoverage.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class DepthOfCoverageNodeModel extends NodeModel {
    
	 protected static final NodeLogger logger = NodeLogger.getLogger(DepthOfCoverageNodeModel.class);
	
	/**
	 * Config Keys
	 */
	
	public static final String CFGKEY_PATH2GATK = "gatk";
	public static final String CFGKEY_PATH2BED = "bed";
	public static final String CFGKEY_INFILE = "infile";
	public static final String CFGKEY_EXTRAFILTERS = "extrafilters";
	public static final String CFGKEY_MEMORY = "memory";
	public static final String CFGKEY_THREADS = "threads";
	public static final String CFGKEY_REFGENOME = "refgenome";
	public static final String CFGKEY_FILESUFFIX = "suffix";
	
	
	
	/**
	 * Node Models
	 */
	private final SettingsModelString m_path2gatk = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_PATH2GATK,"");
	private final SettingsModelOptionalString m_path2bed = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_PATH2BED,"",true);
	private final SettingsModelOptionalString m_extrafilters = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_EXTRAFILTERS,"",false);
	private final SettingsModelIntegerBounded m_memory = new SettingsModelIntegerBounded(DepthOfCoverageNodeModel.CFGKEY_MEMORY, 4, 1,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_threads = new SettingsModelIntegerBounded(DepthOfCoverageNodeModel.CFGKEY_THREADS, 1, 1,Integer.MAX_VALUE);
	private final SettingsModelString m_refgenome = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_REFGENOME,"");
	private final SettingsModelString m_filesuffix = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_FILESUFFIX,"");
	private final SettingsModelString m_infile = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_INFILE,"");

	
    /**
     * Constructor for the node model.
     */
    protected DepthOfCoverageNodeModel() {
    
        super(1, 0);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String outfile = m_infile.getStringValue().replace(".bam","."+m_filesuffix.getStringValue());
    	
    	String cmd = "java -jar -Xmx"+m_memory.getIntValue()+"G " + m_path2gatk.getStringValue();
    	cmd += " -T DepthOfCoverage";
    	cmd += " -R "+m_refgenome.getStringValue();
    	cmd += " -o "+outfile;
    	cmd += " -I "+m_infile.getStringValue();
    	
    	if(!m_path2bed.getStringValue().equals("")){
    		cmd += " -L "+m_path2bed.getStringValue();
    	}
    	
    	cmd += " -nt "+m_threads.getIntValue();
    	
    	if(!m_extrafilters.getStringValue().equals("")){
        	String[] filters = m_extrafilters.getStringValue().split(",");
        	for(String filter : filters){
        		cmd += " --read_filter "+filter;
        	}
    	}

    	DepthOfCoverageNodeModel.logger.info("Running GATK DepthOfCoverage...");
    	DepthOfCoverageNodeModel.logger.info("Log files can be found in "+outfile+".out.log and "+outfile+".err.log");
    	
		Executor.executeCommand(new String[]{cmd}, exec, null, DepthOfCoverageNodeModel.logger, outfile+".out.log", outfile+".err.log", null);
		
    	DepthOfCoverage.processCoverageFile(outfile);
		
        return null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

        // TODO: generated method stub
        return null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_extrafilters.saveSettingsTo(settings);
         m_filesuffix.saveSettingsTo(settings);
         m_infile.saveSettingsTo(settings);
         m_memory.saveSettingsTo(settings);
         m_path2bed.saveSettingsTo(settings);
         m_path2gatk.saveSettingsTo(settings);
         m_refgenome.saveSettingsTo(settings);
         m_threads.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_extrafilters.loadSettingsFrom(settings);
        m_filesuffix.loadSettingsFrom(settings);
        m_infile.loadSettingsFrom(settings);
        m_memory.loadSettingsFrom(settings);
        m_path2bed.loadSettingsFrom(settings);
        m_path2gatk.loadSettingsFrom(settings);
        m_refgenome.loadSettingsFrom(settings);
        m_threads.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_extrafilters.validateSettings(settings);
        m_filesuffix.validateSettings(settings);
        m_infile.validateSettings(settings);
        m_memory.validateSettings(settings);
        m_path2bed.validateSettings(settings);
        m_path2gatk.validateSettings(settings);
        m_refgenome.validateSettings(settings);
        m_threads.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }

}


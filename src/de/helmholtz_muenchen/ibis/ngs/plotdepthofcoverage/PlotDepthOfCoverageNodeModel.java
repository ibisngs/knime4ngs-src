package de.helmholtz_muenchen.ibis.ngs.plotdepthofcoverage;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintWriter;

import org.apache.commons.io.IOUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.FileSearch;

/**
 * This is the model implementation of PlotDepthOfCoverage.
 * 
 *
 * @author 
 */
public class PlotDepthOfCoverageNodeModel extends RNodeModel {
    
	
	/**
	 * Config Keys
	 */
	
	public static final String CFGKEY_INFOLDER = "infolder";
	public static final String CFGKEY_FILESUFFIX = "suffix";
	
	/**
	 * Node Models
	 */
	private final SettingsModelString m_infolder = new SettingsModelString(PlotDepthOfCoverageNodeModel.CFGKEY_INFOLDER,"");
	private final SettingsModelString m_filesuffix = new SettingsModelString(PlotDepthOfCoverageNodeModel.CFGKEY_FILESUFFIX,"");
	
	private static final String SCRIPT_PATH = "ngs" + File.separatorChar + "de" + File.separatorChar + "plotExomeCapture.R";
	
	
    /**
     * Constructor for the node model.
     */
    protected PlotDepthOfCoverageNodeModel() {
    	super(1, 0, SCRIPT_PATH, new String[]{"--infile"}, new String[]{"--output"});
    }

    
    /**
     * {@inheritDoc}
     * @throws Exception 
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{    	
    	    	
    	FileSearch fileSearch = new FileSearch();
    	PrintWriter writer = new PrintWriter(m_infolder.getStringValue()+"/DoCSampleSummary.csv", "UTF-8");	 

    	//try different directory and filename :)
    	fileSearch.searchDirectory(new File(m_infolder.getStringValue()), m_filesuffix.getStringValue());
	    for (String matched : fileSearch.getResult()){
	    	
	    	System.out.println("Found : " + matched);
		
	    	FileInputStream inputStream = new FileInputStream(matched);
	    	try {
	    		String everything = IOUtils.toString(inputStream);
	    		writer.println(matched+"\t"+everything);
	    	} finally {
	    		inputStream.close();
	    	}
	    }
    	writer.close();

    	this.addArgument("--in"   , m_infolder.getStringValue()+"/DoCSampleSummary.csv");
    	
    	
    	super.execute(inData, exec);
    	
    	
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
        	m_filesuffix.saveSettingsTo(settings);
        	m_infolder.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_filesuffix.loadSettingsFrom(settings);
        m_infolder.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_filesuffix.validateSettings(settings);
        m_infolder.validateSettings(settings);
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


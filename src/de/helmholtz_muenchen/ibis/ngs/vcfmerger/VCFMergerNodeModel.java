package de.helmholtz_muenchen.ibis.ngs.vcfmerger;

import java.io.File;
import java.io.IOException;

import net.sf.picard.vcf.MergeVcfs;

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

import de.helmholtz_muenchen.ibis.ngs.vqsr.VQSRNodeModel;

/**
 * This is the model implementation of VCFMerger.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class VCFMergerNodeModel extends NodeModel {
    
    public static final String CFGKEY_INFOLDER					= "infolder";
    public static final String CFGKEY_REGEX						= "regex";
    public static final String CFGKEY_OUTFOLDER					= "outfolder";
    
    private final SettingsModelString m_INFOLDER = new SettingsModelString(VCFMergerNodeModel.CFGKEY_INFOLDER, "");
    private final SettingsModelString m_REGEX = new SettingsModelString(VCFMergerNodeModel.CFGKEY_REGEX, "");
    private final SettingsModelString m_OUTFOLDER = new SettingsModelString(VCFMergerNodeModel.CFGKEY_OUTFOLDER, "");

	
    /**
     * Constructor for the node model.
     */
    protected VCFMergerNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	
    	VCFMerger.mergeVCFs(m_INFOLDER.getStringValue(), m_OUTFOLDER.getStringValue(), m_REGEX.getStringValue(),exec);
    	
        // TODO: Return a BufferedDataTable for each output port 
        return new BufferedDataTable[]{inData[0]};
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
        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_INFOLDER.saveSettingsTo(settings);
         m_OUTFOLDER.saveSettingsTo(settings);
         m_REGEX.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_INFOLDER.loadSettingsFrom(settings);
        m_OUTFOLDER.loadSettingsFrom(settings);
        m_REGEX.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_INFOLDER.validateSettings(settings);
        m_OUTFOLDER.validateSettings(settings);
        m_REGEX.validateSettings(settings);
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


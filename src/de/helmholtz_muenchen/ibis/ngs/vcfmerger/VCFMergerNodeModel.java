package de.helmholtz_muenchen.ibis.ngs.vcfmerger;

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
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;


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
    public static final String CFGKEY_GATK						= "gatk";
	public static final String CFGKEY_REF_GENOME 				= "REFGENOME";
	public static final String CFGKEY_GENOTYPEMERGEOPTION 		= "genotypemergeoption";
    
    private final SettingsModelString m_GATK 				= new SettingsModelString(VCFMergerNodeModel.CFGKEY_GATK, "");
    private final SettingsModelString m_REF_GENOME 			= new SettingsModelString(VCFMergerNodeModel.CFGKEY_REF_GENOME, "");
    private final SettingsModelString m_INFOLDER 			= new SettingsModelString(VCFMergerNodeModel.CFGKEY_INFOLDER, "");
    private final SettingsModelString m_REGEX 				= new SettingsModelString(VCFMergerNodeModel.CFGKEY_REGEX, "");
    private final SettingsModelString m_OUTFOLDER 			= new SettingsModelString(VCFMergerNodeModel.CFGKEY_OUTFOLDER, "");
    private final SettingsModelString m_GENOTYPEMERGEOPTION	= new SettingsModelString(VCFMergerNodeModel.CFGKEY_GENOTYPEMERGEOPTION, "");

	private static final NodeLogger LOGGER = NodeLogger.getLogger(VCFMergerNodeModel.class);
    
	//The Output Col Names
	public static final String OUT_COL1 = "MERGED VARIANTS";
	
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

    	String OUTFILE = VCFMerger.mergeVCFs(m_GATK.getStringValue(),m_REF_GENOME.getStringValue(),m_INFOLDER.getStringValue(), m_OUTFOLDER.getStringValue(), m_REGEX.getStringValue(),m_GENOTYPEMERGEOPTION.getStringValue(),exec,LOGGER);
    	
    	/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(OUTFILE)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();

        return new BufferedDataTable[]{outTable};
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

        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_INFOLDER.saveSettingsTo(settings);
         m_OUTFOLDER.saveSettingsTo(settings);
         m_REGEX.saveSettingsTo(settings);
         m_GATK.saveSettingsTo(settings);
         m_REF_GENOME.saveSettingsTo(settings);
         m_GENOTYPEMERGEOPTION.saveSettingsTo(settings);
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
        m_GATK.loadSettingsFrom(settings);
        m_REF_GENOME.loadSettingsFrom(settings);
        m_GENOTYPEMERGEOPTION.loadSettingsFrom(settings);
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
        m_GATK.validateSettings(settings);
        m_REF_GENOME.validateSettings(settings);
        m_GENOTYPEMERGEOPTION.validateSettings(settings);
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


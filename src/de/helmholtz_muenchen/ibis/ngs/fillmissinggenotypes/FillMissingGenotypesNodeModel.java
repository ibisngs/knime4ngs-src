package de.helmholtz_muenchen.ibis.ngs.fillmissinggenotypes;

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

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;

/**
 * This is the model implementation of FillMissingGenotypes.
 * 
 *
 * @author 
 */
public class FillMissingGenotypesNodeModel extends NodeModel {
    
	
	public static final String CFGKEY_COVERAGEFILEFOLDER 	= "COVERAGEFILEFOLDER";
	public static final String CFGKEY_FILESUFFIX		 	= "FILESUFFIX";
	
    private final SettingsModelString m_COVERAGEFOLDER 		= new SettingsModelString(FillMissingGenotypesNodeModel.CFGKEY_COVERAGEFILEFOLDER, "");
    private final SettingsModelString m_FILESUFFIX 			= new SettingsModelString(FillMissingGenotypesNodeModel.CFGKEY_FILESUFFIX, "");

	//The Output Col Names
	public static final String OUT_COL1 = "FilledGTVARIANTS";
    
    /**
     * Constructor for the node model.
     */
    protected FillMissingGenotypesNodeModel() {
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

//		String PATH_DATA_DIR = "/storageNGS/ngs1/projects/exome/schizo_Rouleau_SCZ/privateData";
//		String FILE_SUFFIX = "filtered_aln_noDups_realigned_recal";	
    	
    	String INFILE  = inData[0].iterator().next().getCell(0).toString();
    	String OUTFILE = FillPindelGenotypes.fillGTs(INFILE, m_COVERAGEFOLDER.getStringValue(), m_FILESUFFIX.getStringValue());
    	
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
    protected void reset() {}

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
         m_COVERAGEFOLDER.saveSettingsTo(settings);
         m_FILESUFFIX.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_COVERAGEFOLDER.loadSettingsFrom(settings);
        m_FILESUFFIX.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_COVERAGEFOLDER.validateSettings(settings);
        m_FILESUFFIX.validateSettings(settings);
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


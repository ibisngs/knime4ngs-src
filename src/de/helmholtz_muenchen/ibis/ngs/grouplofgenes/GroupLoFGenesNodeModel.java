package de.helmholtz_muenchen.ibis.ngs.grouplofgenes;

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
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GroupLoFGenes.
 * 
 *
 * @author tim.jeske
 */
public class GroupLoFGenesNodeModel extends NodeModel {
    	
	static final String CFGKEY_GENESUM = "gene_summary";
	private final SettingsModelString m_genesum = new SettingsModelString(CFGKEY_GENESUM,"-");
	
	static final String CFGKEY_CDSFILE = "cds_file";
	private final SettingsModelString m_cdsfile = new SettingsModelString(CFGKEY_CDSFILE, "-");
	
	static final String CFGKEY_SAMPLES = "samples";
	private final SettingsModelInteger m_samples = new SettingsModelInteger(CFGKEY_SAMPLES,0);
	
	//output col names
	public static final String OUT_COL1 = "Path2NeverObservedGenes";
	public static final String OUT_COL2 = "Path2LoFTolerantGenes";
	public static final String OUT_COL3 = "Path2LoFIntolerantGenes";
	
	public static boolean optionalPort=false;
	
    /**
     * Constructor for the node model.
     */
    protected GroupLoFGenesNodeModel() {
    	super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(1));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	
    	String gene_summary;
    	if(optionalPort) {
    		gene_summary = inData[0].iterator().next().getCell(1).toString();
    	} else {
    		gene_summary = m_genesum.getStringValue();
    	}
    	
    	String cds_file = m_cdsfile.getStringValue();
    	int samples = m_samples.getIntValue();
    	GeneListExtractor gle = new GeneListExtractor(cds_file,gene_summary,samples);
    	String [] outfiles = gle.getGeneLists();
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    					new DataColumnSpec[]{
    	    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    	    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec(),
    	    					new DataColumnSpecCreator(OUT_COL3, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(outfiles[0]),
    			(FileCell) FileCellFactory.create(outfiles[1]),
    			(FileCell) FileCellFactory.create(outfiles[2])};
    	
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

    	try{
			inSpecs[0].getColumnNames();
			optionalPort=true;
			
		}catch(NullPointerException e){
			String warning = CheckUtils.checkSourceFile(m_genesum.getStringValue());
	        if (warning != null) {
	            setWarningMessage(warning);
	        }
		}
    	
    	String warning = CheckUtils.checkSourceFile(m_cdsfile.getStringValue());
        if (warning != null) {
            setWarningMessage(warning);
        }
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL3, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_genesum.saveSettingsTo(settings);
         m_cdsfile.saveSettingsTo(settings);
         m_samples.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_genesum.loadSettingsFrom(settings);
        m_cdsfile.loadSettingsFrom(settings);
        m_samples.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_genesum.validateSettings(settings);
        m_cdsfile.validateSettings(settings);
        m_samples.validateSettings(settings);
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


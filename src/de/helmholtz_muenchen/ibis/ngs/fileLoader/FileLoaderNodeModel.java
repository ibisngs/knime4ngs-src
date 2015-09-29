package de.helmholtz_muenchen.ibis.ngs.fileLoader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

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
 * This is the model implementation of FileLoader.
 * 
 *
 * @author Tim Jeske
 */
public class FileLoaderNodeModel extends NodeModel {
	
	public static final String CFGKEY_INFILE1 = "infile1";
	public static final String CFGKEY_INFILE2 = "infile2";
	
	private final SettingsModelString m_infile1 = new SettingsModelString(FileLoaderNodeModel.CFGKEY_INFILE1,"");
	private final SettingsModelString m_infile2 = new SettingsModelString(FileLoaderNodeModel.CFGKEY_INFILE2,"");
	
	public static final String OUT_COL1 = "Path2File1";
	public static final String OUT_COL2 = "Path2File2";
	
	private static final String [] ENDINGS = {".vcf",".fastq",".bam",".sam"};
	boolean firstOk = false;
	boolean secondOk = false;
	
    /**
     * Constructor for the node model.
     */
    protected FileLoaderNodeModel() {
    	super(0, 1);
    	m_infile2.setEnabled(false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	DataColumnSpec [] colSpec;
    	
    	if(firstOk && secondOk) {
    		colSpec = new DataColumnSpec[]{
    				new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    				new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()};
    	} else {
    		colSpec = new DataColumnSpec[]{
    				new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()};
    	}
    	
    	FileCell [] fileCell;
    	
    	if(firstOk && secondOk) {
    		fileCell = new FileCell[] {
        			(FileCell) FileCellFactory.create(m_infile1.getStringValue()),
        			(FileCell) FileCellFactory.create(m_infile2.getStringValue())};
    	} else {
    		fileCell = new FileCell[] {
        			(FileCell) FileCellFactory.create(m_infile1.getStringValue())};
    	}
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(colSpec));
    	
    	FileCell[] c = fileCell;
    	
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
    	
    	String in1 = m_infile1.getStringValue();
    	String in2 = m_infile2.getStringValue();
    	
    	firstOk = false;
    	secondOk = false;
    	
    	if(Files.notExists(Paths.get(in1)) || !hasAllowedEnding(in1)) {
    		throw new InvalidSettingsException("First input file does not exist or has disallowed ending!");
    	} else {
    		firstOk = true;
    	}
    	
    	if(in2.length()>0) {
    		if(Files.notExists(Paths.get(in2)) || !in2.endsWith(".fastq")) {
    			setWarningMessage("Second input file does not exist or has disallowed ending!");
    		} else {
    			secondOk = true;
    		}
    	}
    	
    	
    	if(firstOk && secondOk) {
    	return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()})};
    	} else {
    		return new DataTableSpec[]{new DataTableSpec(
        			new DataColumnSpec[]{
        					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()})};
    	}
    }
    
    protected boolean hasAllowedEnding(String path) {
    	for(String s: ENDINGS) {
    		if(path.endsWith(s)) {
    			return true;
    		}
    	}
    	return false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_infile1.saveSettingsTo(settings);
         m_infile2.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_infile1.loadSettingsFrom(settings);
        m_infile2.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_infile1.validateSettings(settings);
        m_infile2.validateSettings(settings);
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


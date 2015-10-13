package de.helmholtz_muenchen.ibis.ngs.fileLoader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
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

import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BEDCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.SAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;

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
	
	private static final String [] ENDINGS = {".vcf",".fastq",".bam",".sam", ".bed"};
	private static final DataType [] TYPES = {VCFCell.TYPE, FastQCell.TYPE, BAMCell.TYPE, SAMCell.TYPE, BEDCell.TYPE};
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

    	String in1 = m_infile1.getStringValue();
    	String in2 = m_infile2.getStringValue();
    	
    	DataColumnSpec dcs1 = null;
    	DataColumnSpec dcs2 = null;
    	DataColumnSpec [] specs = null;
    	
    	dcs1 = new DataColumnSpecCreator(OUT_COL1, TYPES[checkEnding(in1)]).createSpec();
		specs = new DataColumnSpec[]{dcs1};

    	if(secondOk) {
			dcs2 = new DataColumnSpecCreator(OUT_COL2, TYPES[checkEnding(in2)]).createSpec();
    		specs = new DataColumnSpec[]{dcs1, dcs2};
    	}
    	
    	FileCell [] fileCell = new FileCell[] {
    			(FileCell) FileCellFactory.create(in1)};
    	
    	if(secondOk) {
    		fileCell = new FileCell[] {
        			(FileCell) FileCellFactory.create(in1),
        			(FileCell) FileCellFactory.create(in2)};
    	}
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(specs));
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
    	
    	DataColumnSpec dcs1 = null;
    	DataColumnSpec dcs2 = null;
    	DataColumnSpec [] specs = null;
    	
    	secondOk = false;
    	
    	//check first input file
    	if(Files.notExists(Paths.get(in1)) || checkEnding(in1)==-1) {
    		throw new InvalidSettingsException("First input file does not exist or has disallowed ending!");
    	}
    	
    	//first input file is ok
    	dcs1 = new DataColumnSpecCreator(OUT_COL1, TYPES[checkEnding(in1)]).createSpec();
    	
    	if(in2.length()>0) {
    		if(Files.notExists(Paths.get(in2)) || !in2.endsWith(".fastq")) {
    			setWarningMessage("Second input file does not exist or has disallowed ending and will be ignored!");
    		} else {
    			secondOk = true;
    			dcs2 = new DataColumnSpecCreator(OUT_COL2, TYPES[checkEnding(in2)]).createSpec();
    		}
    	}
    	
    	if(secondOk) {
    		specs = new DataColumnSpec[]{dcs1, dcs2};
    	} else {
    		specs = new DataColumnSpec[]{dcs1};
    	}
    	
    	return new DataTableSpec[]{new DataTableSpec(specs)};
    }
    
    protected int checkEnding(String path) {
    	for(int i = 0; i < ENDINGS.length; i++) {
    		if(path.endsWith(ENDINGS[i])) {
    			return i;
    		}
    	}
    	return -1;
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


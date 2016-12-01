package de.helmholtz_muenchen.ibis.ngs.Bcl2FastQ;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.lang3.StringUtils;
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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;


/**
 * This is the model implementation of Bcl2FastQ.
 * 
 *
 * @author Kaarin Ahomaa
 */
public class Bcl2FastQNodeModel extends HTExecutorNodeModel {
    
	// keys for SettingsModels
	static final String CFGKEY_TOOL_PATH = "ToolPath";
	final SettingsModelString m_ToolPath = new SettingsModelString(CFGKEY_TOOL_PATH,"");
	
	static final String CFGKEY_INPUT_PATH = "InputPath";
	final SettingsModelString m_InputPath = new SettingsModelString(CFGKEY_INPUT_PATH,"");
	
	static final String CFGKEY_OUTFOLDER = "outfolder";
	final SettingsModelString m_outfolder = new SettingsModelString(CFGKEY_OUTFOLDER,"");

	static final String CFGKEY_ISPAIRED = "IsPaired";
	final SettingsModelBoolean m_IsPaired = new SettingsModelBoolean(CFGKEY_ISPAIRED, true);
	
	static final String CFGKEY_THREADS = "threads";
	final SettingsModelIntegerBounded m_threads = new SettingsModelIntegerBounded(CFGKEY_THREADS,4, 1, Integer.MAX_VALUE);
	
	static final String CFGKEY_OPT_FLAGS = "OptionalFlags";
	final SettingsModelOptionalString m_OptionalFlags = new SettingsModelOptionalString(CFGKEY_OPT_FLAGS,"",false);
	
	
	//static final String CFGKEY_INTEROP = "interop";
	//final SettingsModelString m_interop = new SettingsModelString(CFGKEY_INTEROP, "");
	
	private String OUT_COL1 = "Path2ReadsFile1";
	private String OUT_COL2 = "Path2ReadsFile2";
	
	
	DataColumnSpec out1 = new DataColumnSpecCreator(OUT_COL1, FastQCell.TYPE).createSpec();
	DataColumnSpec out2 = new DataColumnSpecCreator(OUT_COL2, FastQCell.TYPE).createSpec();
	
    /**
     * Constructor for the node model.
     */
    protected Bcl2FastQNodeModel() {
        super(0, 1);
        addSetting(m_ToolPath);
        addSetting(m_InputPath);
    	addSetting(m_outfolder);
    	addSetting(m_IsPaired);
    	addSetting(m_threads);
    	addSetting(m_OptionalFlags);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec)
			throws Exception {

		// tool
		String tool = m_ToolPath.getStringValue();

		// input files
		String infiles = m_InputPath.getStringValue();

		// output files
		String outfiles = m_outfolder.getStringValue();

		// String interop = m_interop.getStringValue();

		int input = m_threads.getIntValue();

		int p = input;
		int r, w, d;

		if (input < 4) {
			r = input;
			w = input;
			d = input;
		} else {
			r = 4;
			w = 4;
			d = (int) Math.round(0.2 * (double) input);
		}

		String cmd = tool;
		cmd += " -R " + infiles;
		cmd += " -o " + outfiles;
		cmd += " --loading-threads " + r;
		cmd += " --demultiplexing-threads " + d;
		cmd += " --processing-threads " + p;
		cmd += " --writing-threads " + w;
//		cmd += " --no-lane-splitting "; 
		cmd += " "+m_OptionalFlags.getStringValue();
		// cmd += " --interop-dir="+ interop;

		String lockFile = outfiles + File.separatorChar + "bcl2fastq" + SuccessfulRunChecker.LOCK_ENDING;
		String stdErrFile = outfiles + File.separatorChar + "bcl2fastq.stderr";
		String stdOutFile = outfiles + File.separatorChar + "bcl2fastq.stdout";

		/** Execute **/
		super.executeCommand(new String[] { StringUtils.join(cmd, " ") }, outfiles, exec, new File(lockFile), stdErrFile,
				stdOutFile);

		File folder = new File(outfiles);
		File[] listOfFiles = folder.listFiles();
		Arrays.sort(listOfFiles);
		File curr;

		ArrayList<String> fastqFiles = new ArrayList<>();

		for (int i = 0; i < listOfFiles.length; i++) {
			curr = listOfFiles[i];
			if (curr.isFile() && curr.getName().endsWith(".fastq.gz")) {
				fastqFiles.add(curr.getAbsolutePath());
			}
		}
		


		// Create Output Table
		BufferedDataContainer cont = exec.createDataContainer(new DataTableSpec(new DataColumnSpec[] { out1 }));

		FileCell[] c;

		if (m_IsPaired.getBooleanValue()) {
			String warning;
			if (fastqFiles.size() % 2 == 1) {
				warning = "Odd number of fastq files. Check wheather you have paired end reads in your output directory.";
				setWarningMessage(warning);
			}
			cont = exec.createDataContainer(new DataTableSpec(new DataColumnSpec[] { out1, out2 }));

			for (int i = 0; i < fastqFiles.size(); i = i + 2) {
				c = new FileCell[] { (FileCell) FileCellFactory.create(fastqFiles.get(i)),
						(FileCell) FileCellFactory.create(fastqFiles.get(i + 1)) };
				cont.addRowToTable(new DefaultRow("Row" + i / 2, c));
			}

		} else {
			for (int i = 0; i < fastqFiles.size(); i++) {
				c = new FileCell[] { (FileCell) FileCellFactory.create(fastqFiles.get(i)) };
				cont.addRowToTable(new DefaultRow("Row" + i, c));

			}
		}

		cont.close();
		BufferedDataTable outTable = cont.getTable();

		return new BufferedDataTable[] { outTable };
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

        if(m_IsPaired.getBooleanValue()) {
        	return new DataTableSpec[]{new DataTableSpec(
        			new DataColumnSpec[]{out1,out2})};
        }
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{out1})};
    }
   
    
//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    protected void saveSettingsTo(final NodeSettingsWO settings) {
//         // TODO: generated method stub
//    	/** added for HTE **/
//    	super.saveSettingsTo(settings);
//    	
//    	m_ToolPath.saveSettingsTo(settings);
//    	m_InputPath.saveSettingsTo(settings);
//    	m_outfolder.saveSettingsTo(settings);
//    	m_IsPaired.saveSettingsTo(settings);
//    	m_threads.saveSettingsTo(settings);
//    	//m_interop.saveSettingsTo(settings);
//    }
//
//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
//            throws InvalidSettingsException {
//        // TODO: generated method stub
//    	/** added for HTE **/
//    	super.loadValidatedSettingsFrom(settings);
//    	
//    	m_ToolPath.loadSettingsFrom(settings);
//    	m_InputPath.loadSettingsFrom(settings);
//    	m_outfolder.loadSettingsFrom(settings);
//    	m_IsPaired.loadSettingsFrom(settings);
//    	m_threads.loadSettingsFrom(settings);
//    	//m_interop.loadSettingsFrom(settings);
//    	
//    }
//
//    /**
//     * {@inheritDoc}
//     */
//    @Override
//    protected void validateSettings(final NodeSettingsRO settings)
//            throws InvalidSettingsException {
//        // TODO: generated method stub
//    	/** added for HTE **/
//    	super.validateSettings(settings);
//    	
//    	m_ToolPath.validateSettings(settings);
//    	m_InputPath.validateSettings(settings);
//    	m_outfolder.validateSettings(settings);
//    	m_IsPaired.validateSettings(settings);
//    	m_threads.validateSettings(settings);
//    	//m_interop.validateSettings(settings);
//    	
//    }
    
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


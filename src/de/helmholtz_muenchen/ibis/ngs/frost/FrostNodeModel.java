package de.helmholtz_muenchen.ibis.ngs.frost;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.RowKey;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.IntCell;
import org.knime.core.data.def.StringCell;
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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.frost.FrostRunner;
import de.helmholtz_muenchen.ibis.utils.ngs.frost.Trio_Simulator;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of Frost.
 * 
 *
 * @author Syeda Tanzeem Haque
 */
public class FrostNodeModel extends NodeModel {
	// the logger instance
    private static final NodeLogger LOGGER = NodeLogger.getLogger(FrostNodeModel.class); //not used yet
    private int seed = -11;
    
	public static final String INTERNAL_PATH = "/home/ibis/tanzeem.haque/Documents/Frost_outputs/";

	public static final String CFGKEY_FASTA = "fasta";
	public static final String CFGKEY_MUT_RATE = "mutation";
	public static final String CFGKEY_RECOMB_NUM = "recombination";
	public static final String CFGKEY_SEED = "seed";
	
	public static final String CFGKEY_USE_MUT_RATE = "use_mutation";
	public static final String CFGKEY_USE_RECOMB_NUM = "use_recombination";
	public static final String CFGKEY_USE_SEED = "use_seed";



    /** initial default count value. */
    static final double DEFAULT_MUTATION_RATE = 2.36;
    static final int DEFAULT_RECOMNATION = 1000;
    static final int DEFAULT_SEED = 999;

    private final SettingsModelString m_FASTA = new SettingsModelString(FrostNodeModel.CFGKEY_FASTA,""); // path to fasta file
    private final SettingsModelDoubleBounded m_MUT = new SettingsModelDoubleBounded(FrostNodeModel.CFGKEY_MUT_RATE,FrostNodeModel.DEFAULT_MUTATION_RATE,1.0, 3.0);
    private final SettingsModelIntegerBounded m_REC = new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_RECOMB_NUM,FrostNodeModel.DEFAULT_RECOMNATION,0, Integer.MAX_VALUE);
    private final SettingsModelIntegerBounded m_SEED = new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_SEED,FrostNodeModel.DEFAULT_SEED,Integer.MIN_VALUE, Integer.MAX_VALUE);
	
    /**
	 * copy the booleans into dialog class to add component and check box config
	 */
    
    private final SettingsModelBoolean m_use_MUT = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_MUT_RATE, false);
	private final SettingsModelBoolean m_use_REC = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_RECOMB_NUM, false);
	private final SettingsModelBoolean m_use_SEED = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_SEED, false);

	/**The Output Col Names */
	public static final String OUT_COL = "fasta_output_";
	public static final String OUT_COL_LAST = "record_files";
    /**
     * Constructor for the node model.
     */
    protected FrostNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(0, 1);
        m_MUT.setEnabled(false);
        m_REC.setEnabled(false);
        m_SEED.setEnabled(false);
    }

    
    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

//		Runner.run(input, mutRate, recombination, seed);
    	/**
    	 * Parameters
    	 */
    	String input = m_FASTA.getStringValue();
        double mutRate = 1.;
        int rec_num = -1;
        /**
         * the seed is set here and it runs the FrostRunner.main(args)
         */
        executeInput(input, mutRate, rec_num); 
    	/**
    	 * OUTPUT
    	 */
    	 // the data table spec of the single output table, 
        // the table will have three columns:
        
        int col_num = 1 + FrostRunner.ID_List.length;
        DataColumnSpec[] allColSpecs = new DataColumnSpec[col_num];
        /**
         * create the last column which will have 6 rows: record files
         */
        allColSpecs[allColSpecs.length-1] = new DataColumnSpecCreator(OUT_COL_LAST, FileCell.TYPE).createSpec();
        /**
         * create the column(s) for (multi)fasta output of a trio with 6 rows
         */
    	for (int i = 0; i < FrostRunner.ID_List.length; i++) {
            allColSpecs[i] = new DataColumnSpecCreator("Col. " + (i+1) + ": " + OUT_COL + FrostRunner.ID_List[i], FileCell.TYPE).createSpec();

    	}
        DataTableSpec outputSpec = new DataTableSpec(allColSpecs);
        BufferedDataContainer container = exec.createDataContainer(outputSpec);
        /**
         * creating the 6 rows
         * i = 0,1,2,3,4,5
         */
        FileCell[] cells = new FileCell[col_num];
        for (int i = 0; i < 6; i++) {
            RowKey key = new RowKey("Row " + (i+1));
            // the cells of the current row, the types of the cells must match
            //number of cells = col_num
            // the column spec (see above)
			cells[cells.length - 1] = (FileCell) FileCellFactory.create(recordFiles()[i]);
			
			for (int j = 0; j < FrostRunner.ID_List.length; j++) {
				switch (i) {
				case 0:
					cells[j] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_PATH + FrostRunner.ID_List[j] + "_F_0.fa");
					break;
				case 1:
					cells[j] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_PATH + FrostRunner.ID_List[j] + "_F_1.fa");
					break;
				case 2:
					cells[j] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_PATH + FrostRunner.ID_List[j] + "_M_0.fa");
					break;
				case 3:
					cells[j] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_PATH + FrostRunner.ID_List[j] + "_M_1.fa");
					break;
				case 4:
					cells[j] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_PATH + FrostRunner.ID_List[j] + "_C_0.fa");
					break;
				case 5:
					cells[j] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_PATH + FrostRunner.ID_List[j] + "_C_1.fa");
					break;
				}

			}
            DataRow row = new DefaultRow(key, cells);
            container.addRowToTable(row);
            
            // check if the execution monitor was canceled
            exec.checkCanceled();        

        }
        // once we are done, we close the container and return its table
        container.close();
        BufferedDataTable out = container.getTable();
        return new BufferedDataTable[]{out};
    }
    
    private String[] recordFiles() {
    	
    	String[] files = new String[6];
    	files[0] = FrostNodeModel.INTERNAL_PATH + "ids.txt";
    	files[1] = FrostNodeModel.INTERNAL_PATH + "deNovo_" + this.seed + ".txt";
    	files[2] = FrostNodeModel.INTERNAL_PATH + "recombination_" + this.seed + ".txt";
    	files[3] = FrostNodeModel.INTERNAL_PATH + "parents_run_" + this.seed + ".txt";
    	files[4] = FrostNodeModel.INTERNAL_PATH + "child_run_" + this.seed + ".txt";
    	files[5] = FrostNodeModel.INTERNAL_PATH + "recombined_seq_" + this.seed + ".txt";
    	return files;

    }
    
    private void executeInput(String input, double mutRate, int rec_num) throws InterruptedException {
        
        ArrayList<String> command = new ArrayList<String>();
    	command.add("-i");
    	command.add(input);

    	if(m_use_MUT.getBooleanValue()) {
    		mutRate = m_MUT.getDoubleValue();
    		command.add("-m");
        	command.add(mutRate+"");
    	}
    	if(m_use_REC.getBooleanValue()) {
    		rec_num = m_REC.getIntValue();
    		command.add("-r");
        	command.add(rec_num+"");
    	}
    	if(m_use_SEED.getBooleanValue()){
    		this.seed = m_SEED.getIntValue();
    		command.add("-s");
        	command.add(seed+"");
    	}
    	
    	String[] args = new String[command.size()];
    	args = command.toArray(args);
//    	String [] args = (String[]) command.toArray();
    	FrostRunner.main(args);

    	
        /**
         * Optional tags -m -r -s
         */
    	
//		String command = "-i "+ input +" -m "+ mutRate + " -r " + rec_num + " -s " + seed; 
    
    	/**Execute**/
//    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
			
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

//    	throw new InvalidSettingsException(s);
    	
        // TODO: generated method stub
        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         // TODO: generated method stub
    	m_FASTA.saveSettingsTo(settings);
    	m_MUT.saveSettingsTo(settings);
    	m_REC.saveSettingsTo(settings);
    	m_SEED.saveSettingsTo(settings);
    	m_use_MUT.saveSettingsTo(settings);
    	m_use_REC.saveSettingsTo(settings);
    	m_use_SEED.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        // TODO: generated method stub
    	m_FASTA.loadSettingsFrom(settings);
    	m_MUT.loadSettingsFrom(settings);
    	m_REC.loadSettingsFrom(settings);
    	m_SEED.loadSettingsFrom(settings);
    	m_use_MUT.loadSettingsFrom(settings);
    	m_use_REC.loadSettingsFrom(settings);
    	m_use_SEED.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        // TODO: generated method stub
    	m_FASTA.validateSettings(settings);
    	m_MUT.validateSettings(settings);
    	m_REC.validateSettings(settings);
    	m_SEED.validateSettings(settings);
    	m_use_MUT.validateSettings(settings);
    	m_use_REC.validateSettings(settings);
    	m_use_SEED.validateSettings(settings);
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


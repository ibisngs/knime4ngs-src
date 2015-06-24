package de.helmholtz_muenchen.ibis.ngs.frost;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.RowKey;
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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.frost.FrostRunner;

/**
 * This is the model implementation of Frost.
 * 
 *
 * @author Syeda Tanzeem Haque
 */
public class FrostNodeModel extends NodeModel {
	// the logger instance
//    private static final NodeLogger LOGGER = NodeLogger.getLogger(FrostNodeModel.class); //not used yet
    private static int seed;
	public static String exome = "";

    
	public static final String INTERNAL_OUTPUT_PATH = "/storageNGS/scratch/Sequenciator/Frost/";
//	public static final String INTERNAL_OUTPUT_PATH = "/home/ibis/tanzeem.haque/Documents/Frost_outputs/";


	public static final String CFGKEY_FASTA = "fasta";
	public static final String CFGKEY_MUT_RATE = "mutation";
	public static final String CFGKEY_RECOMB_NUM = "recombination";
	public static final String CFGKEY_GENERATION = "generation";
	public static final String CFGKEY_SEED = "seed";
	public static final String CFGKEY_VARY="parameter";
	public static final String CFGKEY_BED_FILE="bedfile";


	
	public static final String CFGKEY_USE_MUT_RATE = "use_mutation";
	public static final String CFGKEY_USE_RECOMB_NUM = "use_recombination";
	public static final String CFGKEY_USE_GENERATION = "use_generation";
	public static final String CFGKEY_USE_SEED = "use_seed";
	public static final String CFGKEY_USE_BED_FILE = "use_bedfile";




    /** initial default count value. */
    static final double DEFAULT_MUTATION_RATE = 2.36;
    static final int DEFAULT_RECOMNATION = 1000;
    static final int DEFAULT_GENERATION = 5300;
    static final int DEFAULT_SEED = 999;
    public final static String DEFAULT_MAPFILE = "/storageNGS/scratch/Sequenciator/secondary_files/UCSC_CodingExons.bed";

    private final SettingsModelString m_FASTA = new SettingsModelString(FrostNodeModel.CFGKEY_FASTA,""); // path to fasta file
    private final SettingsModelDoubleBounded m_MUT = new SettingsModelDoubleBounded(FrostNodeModel.CFGKEY_MUT_RATE,FrostNodeModel.DEFAULT_MUTATION_RATE,1.0, 3.0);
    private final SettingsModelIntegerBounded m_REC = new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_RECOMB_NUM,FrostNodeModel.DEFAULT_RECOMNATION,0, Integer.MAX_VALUE);
    private final SettingsModelIntegerBounded m_GEN = new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_GENERATION,FrostNodeModel.DEFAULT_GENERATION,0, Integer.MAX_VALUE);
    private final SettingsModelIntegerBounded m_SEED = new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_SEED,FrostNodeModel.DEFAULT_SEED,Integer.MIN_VALUE, Integer.MAX_VALUE);
    private final SettingsModelString m_VARY = new SettingsModelString(FrostNodeModel.CFGKEY_VARY, "--");
    private final SettingsModelString m_BED_FILE = new SettingsModelString(FrostNodeModel.CFGKEY_BED_FILE,""); // path to fasta file

    /**
	 * copy the booleans into dialog class to add component and check box config
	 */
    
    private final SettingsModelBoolean m_use_MUT = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_MUT_RATE, false);
	private final SettingsModelBoolean m_use_REC = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_RECOMB_NUM, false);
	private final SettingsModelBoolean m_use_GEN = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_GENERATION, false);
	private final SettingsModelBoolean m_use_SEED = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_SEED, false);
	private final SettingsModelBoolean m_use_BED_FILE = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_BED_FILE, false);
	

	/**The Output Row Names */
	public static final String OUT_COL = "fasta_output_";
	/* public static final String OUT_ROW_LAST = "record_files_frost";*/
    /**
     * Constructor for the node model.
     */
    protected FrostNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(0, 1);
//        m_use_MUT.setEnabled(false);
        m_MUT.setEnabled(false);
       
//        m_use_REC.setEnabled(false);
        m_REC.setEnabled(false);
        
//        m_use_GEN.setEnabled(false);
        m_GEN.setEnabled(false);
        
//        m_use_SEED.setEnabled(false);
        m_SEED.setEnabled(false);
        
//        m_use_BED_FILE.setEnabled(false);
        m_BED_FILE.setEnabled(false);

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
        int gen_num = -1;
//        String param = m_VARY.getStringValue();
        /**
         * the seed is set here and it runs the FrostRunner.main(args)
         */
        executeInput(input, mutRate, rec_num, gen_num); 
    	/**
    	 * OUTPUT
    	 */
    	 // the data table spec of the single output table, 
        // the table will have six columns:
        
        
        int col_num = 6; // 6 fastas
        DataColumnSpec[] allColSpecs = new DataColumnSpec[col_num];
        
        /**
         * create the column(s) for (multi)fasta output of a trio
         */ 
		allColSpecs[0] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "M0.fa", FileCell.TYPE).createSpec();
		allColSpecs[1] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "M1.fa", FileCell.TYPE).createSpec();
		allColSpecs[2] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "F0.fa", FileCell.TYPE).createSpec();
		allColSpecs[3] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "F1.fa", FileCell.TYPE).createSpec();
		allColSpecs[4] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "C0.fa", FileCell.TYPE).createSpec();
		allColSpecs[5] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "C1.fa", FileCell.TYPE).createSpec();
		
		
        DataTableSpec outputSpec = new DataTableSpec(allColSpecs);
        BufferedDataContainer container = exec.createDataContainer(outputSpec);       
        /**
         * creating the rows
         * i = 0,1, ... #IDS
         */  
		FileCell[] cells = new FileCell[6];
//		System.out.println("ID LIST now: " + FrostRunner.id_list.size());
		for (int j = 0; j < FrostRunner.id_list.size(); j++) {
			RowKey key = new RowKey("Row: " + FrostRunner.id_list.get(j));
//			System.out.println(key.toString());
			
			for (int i = 0; i < 6 ; i++) {
			// the cells of the current row, the types of the cells must match
			// number of cells = col_num
			// the column spec (see above)

				switch (i) {
					case 0:
						cells[i] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_OUTPUT_PATH + FrostRunner.id_list.get(j) + "_M0.fa");
						break;
					case 1:
						cells[i] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_OUTPUT_PATH + FrostRunner.id_list.get(j) + "_M1.fa");
						break;
					case 2:
						cells[i] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_OUTPUT_PATH + FrostRunner.id_list.get(j) + "_F0.fa");
						break;
					case 3:
						cells[i] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_OUTPUT_PATH + FrostRunner.id_list.get(j) + "_F1.fa");
						break;
					case 4:
						cells[i] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_OUTPUT_PATH + FrostRunner.id_list.get(j) + "_C0.fa");
						break;
					case 5:
						cells[i] = (FileCell) FileCellFactory.create(FrostNodeModel.INTERNAL_OUTPUT_PATH + FrostRunner.id_list.get(j) + "_C1.fa");
						break;
				}
					
			}
			DataRow row = new DefaultRow(key, cells);
            container.addRowToTable(row);

		}
		FrostRunner.id_list.clear();
		/**
         * Record files as flow variable
         */
        /**
         * remember to change the names of parents_run: recordfile[0] and child_run: recordfile[1] to .txt
         */
		String[] fs = FrostNodeModel.recordFiles();
		
		fs[0] = fs[0].replace(".tmp", ".txt");
		fs[1] = fs[1].replace(".tmp", ".txt");
        for(int i = 0; i < fs.length; i++) {
            pushFlowVariableString("record file " +(i+1), fs[i]);
		}
        
//		System.out.println("ID LIST after : " + FrostRunner.id_list.size());

		
		/**
		 * RowKey key = new RowKey("Row: " + FrostNodeModel.OUT_ROW_LAST);
		 * for(int i = 0; i < 6; i++) {
			cells[i] = (FileCell) FileCellFactory.create(recordFiles()[i]);
			System.out.println(cells[i].toString());
			}
		 * DataRow row = new DefaultRow(key, cells);
		 * container.addRowToTable(row);
		 */
	    // check if the execution monitor was canceled
        exec.checkCanceled();
		// once we are done, we close the container and return its table
        container.close();
        BufferedDataTable out = container.getTable();
        return new BufferedDataTable[]{out};
        
        
    }
//    private void createDirectory() {
//		// TODO Auto-generated method stub
//    	File theDir = new File("new folder");
//
//    	// if the directory does not exist, create it
//    	if (!theDir.exists()) {
//    	    System.out.println("creating directory: " + directoryName);
//    	    boolean result = false;
//
//    	    try{
//    	        theDir.mkdir();
//    	        result = true;
//    	    } 
//    	    catch(SecurityException se){
//    	        //handle it
//    	    }        
//    	    if(result) {    
//    	        System.out.println("DIR created");  
//    	    }
//    	}
//		
//	}

    public static String[] recordFiles() {
    	

    	String[] files = new String[7];
		files[0] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "parents_run_" + FrostNodeModel.seed + ".tmp";
    	files[1] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "child_run_" + FrostNodeModel.seed + ".tmp";
    	files[2] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "recombination_" + FrostNodeModel.seed + ".txt";
    	files[3] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "recombined_seq_" + FrostNodeModel.seed + ".txt";
    	files[4] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "trio_phased_" + FrostNodeModel.seed + ".vcf";
    	files[5] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "trio_unphased_" + FrostNodeModel.seed + ".vcf";
    	files[6] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "MyLogFile_" + FrostNodeModel.seed + ".log";

    	return files;

    }
        
    private void executeInput(String input, double mutRate, int rec_num, int gen_num) throws InterruptedException, IOException {
        
    	String param = "";
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
    	if(m_use_GEN.getBooleanValue()){
    		gen_num = m_GEN.getIntValue();
    		command.add("-g");
        	command.add(gen_num+"");
    	}
    	if(m_use_SEED.getBooleanValue()){
    		FrostNodeModel.seed = m_SEED.getIntValue();
    		command.add("-s");
        	command.add(seed+"");
    	}
    	if(m_use_BED_FILE.getBooleanValue()){
			command.add("--bed");
			command.add(m_BED_FILE.getStringValue());
    	}
    	if(m_VARY.getStringValue().equals("Mutation")){
    		param = "--mut";
    	}
    	if(m_VARY.getStringValue().equals("Crossover")){
    		param = "--reco";
    	}
    	if(m_VARY.getStringValue().equals("Denovo")){
    		param = "--denovo";
    	}
    	command.add(param);
    	
//    	for (String s: command)
//    		System.out.println(s);
    	
    	String[] args = new String[command.size()];
    	args = command.toArray(args);
//    	String [] args = (String[]) command.toArray();
    	FrostRunner.main(args);

    	
        /**
         * Optional tags -m -r -s
         */
    	
//		String command = "-i "+ input +" -m "+ mutRate + " -r " + rec_num + " -s " + seed; 
    
			
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

    	int col_num = 6; // 6 fastas
        DataColumnSpec[] allColSpecs = new DataColumnSpec[col_num];
        /**
         * create the column(s) for (multi)fasta output of a trio
         */ 
		allColSpecs[0] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "M0.fa", FileCell.TYPE).createSpec();
		allColSpecs[1] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "M1.fa", FileCell.TYPE).createSpec();
		allColSpecs[2] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "F0.fa", FileCell.TYPE).createSpec();
		allColSpecs[3] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "F1.fa", FileCell.TYPE).createSpec();
		allColSpecs[4] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "C0.fa", FileCell.TYPE).createSpec();
		allColSpecs[5] = new DataColumnSpecCreator(FrostNodeModel.OUT_COL + "C1.fa", FileCell.TYPE).createSpec();
		
		
        DataTableSpec outputSpec = new DataTableSpec(allColSpecs);
    	
        // TODO: generated method stub
        return new DataTableSpec[]{outputSpec};
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
    	m_GEN.saveSettingsTo(settings);
    	m_SEED.saveSettingsTo(settings);
    	m_BED_FILE.saveSettingsTo(settings);
    	m_VARY.saveSettingsTo(settings);
    	m_use_MUT.saveSettingsTo(settings);
    	m_use_REC.saveSettingsTo(settings);
    	m_use_GEN.saveSettingsTo(settings);
    	m_use_SEED.saveSettingsTo(settings);
    	m_use_BED_FILE.saveSettingsTo(settings);
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
    	m_GEN.loadSettingsFrom(settings);
    	m_SEED.loadSettingsFrom(settings);
    	m_BED_FILE.loadSettingsFrom(settings);
    	m_VARY.loadSettingsFrom(settings);
    	m_use_MUT.loadSettingsFrom(settings);
    	m_use_REC.loadSettingsFrom(settings);
    	m_use_GEN.loadSettingsFrom(settings);
    	m_use_SEED.loadSettingsFrom(settings);
    	m_use_BED_FILE.loadSettingsFrom(settings);
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
    	m_GEN.validateSettings(settings);
    	m_SEED.validateSettings(settings);
    	m_BED_FILE.validateSettings(settings);
    	m_VARY.validateSettings(settings);
    	m_use_MUT.validateSettings(settings);
    	m_use_REC.validateSettings(settings);
    	m_use_GEN.validateSettings(settings);
    	m_use_SEED.validateSettings(settings);
    	m_use_BED_FILE.validateSettings(settings);
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


package de.helmholtz_muenchen.ibis.ngs.art;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.RowKey;
import org.knime.core.data.container.CloseableRowIterator;
import org.knime.core.data.def.DefaultRow;

import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;
import de.helmholtz_muenchen.ibis.utils.ngs.frost.FrostRunner;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;
import de.helmholtz_muenchen.ibis.utils.threads.UnsuccessfulExecutionException;
import de.helmholtz_muenchen.ibis.ngs.frostslilhelper.FrostsLilHelperNodeModel;
import de.helmholtz_muenchen.ibis.ngs.frost.FrostNodeModel;


/**
 * This is the model implementation of Art.
 * 
 *
 * @author Syeda Tanzeem Haque
 */
public class ArtNodeModel extends NodeModel {
    
	static boolean optionalPort = false;
	private static final NodeLogger LOGGER = NodeLogger.getLogger(ArtNodeModel.class); //not used yet
    

	public static final String ART_PATH = "/home/ibis/tanzeem.haque/Documents/3rd_party_tools/art/art_bin_VanillaIceCream/art_illumina";
//	public static final String INTERNAL_OUTPUT_PATH = "/storageNGS/scratch/Sequenciator/chr21_art/"; 
	public static final String INTERNAL_OUTPUT_PATH = "/storageNGS/scratch/Sequenciator/hg19_art/"; 

//	public static final String INTERNAL_OUTPUT_PATH = "/storageNGS/scratch/Sequenciator/Art/"; 
//	public static final String INTERNAL_OUTPUT_PATH = "/home/ibis/tanzeem.haque/Documents/Art_outputs/"; 


	/**
	 * is activated iff optionalPort is set to false
	 */
	public static final String CFGKEY_ID_PATH = "id";
	
	/**
//	 *-l the length of reads to be simulated
	 */
	public static final String CFGKEY_LENGTH = "length";
	/**
	 * -m the mean size of DNA/RNA fragments for paired-end simulations
	 */
	public static final String CFGKEY_MEAN_SIZE = "mean_size";
	/**
	 * -f the fold of read coverage to be simulated or number of reads/read pairs generated for each amplicon
	 */
	public static final String CFGKEY_FOLD = "fold";
		/**
	 * -s the standard deviation of DNA/RNA fragment size for paired-end simulations.
	 */
	public static final String CFGKEY_SD = "standard_dev";
	/**
	 * Will be used only as checkbox: optional parameters
	 */
	public static final String CFGKEY_USE_NO_MASK = "use_nf_0";
	public static final String CFGKEY_USE_MATE_PAIR = "use_mp";
	public static final String CFGKEY_USE_ERROR_FREE = "use_ef";
	public static final String CFGKEY_USE_SEPERATE_PROFILE = "use_sp";
	public static final String CFGKEY_USE_QUIET = "use_q";

    /** initial default count value. */
    static final int DEFAULT_LENGTH = 50;
    static final int DEFAULT_MEAN_SIZE = 200;
    static final int DEFAULT_FOLD = 10;
    static final int DEFAULT_SD = 10;

    private final SettingsModelString m_ID_PATH = new SettingsModelString(ArtNodeModel.CFGKEY_ID_PATH,""); // path to fasta files
    /**
     * User defined parameters
     */
    private final SettingsModelIntegerBounded m_LENGTH = new SettingsModelIntegerBounded(ArtNodeModel.CFGKEY_LENGTH, ArtNodeModel.DEFAULT_LENGTH,1, 250);
    private final SettingsModelIntegerBounded m_MEAN_SIZE = new SettingsModelIntegerBounded(ArtNodeModel.CFGKEY_MEAN_SIZE, ArtNodeModel.DEFAULT_MEAN_SIZE,1, Integer.MAX_VALUE);
    private final SettingsModelIntegerBounded m_FOLD = new SettingsModelIntegerBounded(ArtNodeModel.CFGKEY_FOLD, ArtNodeModel.DEFAULT_FOLD,1, 100);
    private final SettingsModelIntegerBounded m_SD = new SettingsModelIntegerBounded(ArtNodeModel.CFGKEY_SD, ArtNodeModel.DEFAULT_LENGTH,1, 100);

    /**
	 * check box config
	 */   
    private final SettingsModelBoolean m_use_NO_MASK = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_NO_MASK, false); //if true, use -nf 0
	private final SettingsModelBoolean m_use_MATE_PAIR = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_MATE_PAIR, false); //if true, use -mp: if -m >2000 -mp on
	private final SettingsModelBoolean m_use_ERROR_FREE = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_ERROR_FREE, false); //if true, use -ef
	private final SettingsModelBoolean m_use_SEPERATE_PROFILE = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_SEPERATE_PROFILE, false); //if true, use -sp
	private final SettingsModelBoolean m_use_QUIET = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_QUIET, false); //if false, do not use -q


	/**The Output Col Names */
	public static final String OUT_COL = "fastq_output_";
	public static ArrayList<String> IDS = new ArrayList<>();
	private ArrayList<Art_object> art_arrList = new ArrayList<Art_object>(); //May be there are more than one chromosome

	

    /**
     * Constructor for the node model.
     */
    protected ArtNodeModel() {
    
        super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(1));
        /**
         * the user defined param-values are enabled
         */
        m_LENGTH.setEnabled(true);
        m_MEAN_SIZE.setEnabled(true);
        m_FOLD.setEnabled(true);
        m_SD.setEnabled(true);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws CanceledExecutionException, InterruptedException,
			ExecutionException, UnsuccessfulExecutionException, IOException {

        // TODO do something here
//        LOGGER.info("Node Model Stub... this is not yet implemented !");
    	/**
    	 * Handling the inputs
    	 */	
    	if(optionalPort){
    		CloseableRowIterator it = inData[0].iterator();
//    		String record = inData[0].getSpec().getName();
//    		System.out.println(record);
    		while (it.hasNext()) {
    			DataRow row = it.next();
    			String id = row.getKey().toString().substring("Row: ".length());
//    			String id = fasta_info[0];
    			/*if (id.startsWith("record_files")) {
    				System.err.println(id);
    				break;
    			}
    			System.out.println("ART KNIME ID: " + id);*/
        		/**
        		 * ArtNodeModel.IDS = new ArrayList<String>();
        		 * ArtNodeModel.IDS.add(id);
        		 */
    			
    			String[] trio = {
    					row.getCell(0).toString(), //M0
    					row.getCell(1).toString(), //M1
    					row.getCell(2).toString(), //F0
    					row.getCell(3).toString(), //F1
    					row.getCell(4).toString(), //C0
    					row.getCell(5).toString() // C1
    					};
//    			for (String s: trio)
//    				System.out.println("inData 0: "+s);
    			/**
    			 * Build a Knime datatyp for the trios (TODO)
    			 */
        		art_arrList = new ArrayList<Art_object>();
    			art_arrList.add(new Art_object(id, trio));
    			if (art_arrList.size() > 0)
    				return output(exec);

    		}

    	}
    	else{
    		String path2ids = m_ID_PATH.getStringValue();
    		/**
    		 * ArtNodeModel.IDS = new ArrayList<String>();
    		 * ArtNodeModel.IDS = getIds(path2ids);
    		 */
    		art_arrList = getUserInput(path2ids, getIds(path2ids));
    	}   
    	return output(exec); 
    	
    }

	private BufferedDataTable[] output(final ExecutionContext exec)
			throws CanceledExecutionException, InterruptedException, ExecutionException, UnsuccessfulExecutionException, IOException {
		
		/**
    	 * here we go...
    	 */
		long startTime = System.currentTimeMillis();

    	execute_help(exec, art_arrList);
    	
    	long endTime   = System.currentTimeMillis();
		NumberFormat formatter = new DecimalFormat("#0.00000");
		System.out.println("Execution time is (art) " + formatter.format((endTime - startTime) / 1000d) + " seconds");
		
//    	System.out.println("REAL SIZE: " + art_arrList.size());

		int col_num = 2; //ArtNodeModel.IDS.size();
    	DataColumnSpec[] allColSpecs = new DataColumnSpec[col_num];/**
    	 * create the column(s) for (multi)fastq output of a trio with 6 rows
    	 */
    	/**
    	for (int i = 0; i < col_num; i++) {
    	    allColSpecs[i] = new DataColumnSpecCreator("Col. " + (i+1) + ": " + OUT_COL + ArtNodeModel.IDS.get(i), FileCell.TYPE).createSpec();

    	}**/
    	allColSpecs[0] = new DataColumnSpecCreator("Col. 1: " + OUT_COL + "fwd", FileCell.TYPE).createSpec();
    	allColSpecs[1] = new DataColumnSpecCreator("Col. 2: " + OUT_COL + "rev", FileCell.TYPE).createSpec();

    	DataTableSpec outputSpec = new DataTableSpec(allColSpecs);
    	BufferedDataContainer container = exec.createDataContainer(outputSpec);
    	// let's add m_count rows to it
    	/**
    	 * creating the 6 rows
    	 * i = 0,1,2,3,4,5
    	 */
    	FileCell[] cells = new FileCell[col_num];

    	int row_idx = 0;
    	for (int i = 0; i < art_arrList.size(); i++) {
    		String chr = art_arrList.get(i).getId();
//    		System.out.println("IDs for output: " + chr);
    	    RowKey key = new RowKey("Row row row your boat");

    		for (int j = 0; j < 3; j++) {
    	        switch (j) {
    			case 0:
    	            key = new RowKey("Row " + (row_idx++) + ". " + chr + "_M");
//    	            System.out.println("key 0: "+ key.toString());
    				cells[0] = (FileCell) FileCellFactory.create(ArtNodeModel.INTERNAL_OUTPUT_PATH + chr + "_M_fwd.fq");
    				cells[1] = (FileCell) FileCellFactory.create(ArtNodeModel.INTERNAL_OUTPUT_PATH + chr + "_M_rev.fq");
    				DataRow row = new DefaultRow(key, cells);
    	    	    container.addRowToTable(row);
    				break;
    			case 1:
    	            key = new RowKey("Row " + (row_idx++) + ". " + chr + "_F");
//    	            System.out.println("key 1: "+ key.toString());
    				cells[0] = (FileCell) FileCellFactory.create(ArtNodeModel.INTERNAL_OUTPUT_PATH + chr + "_F_fwd.fq");
    				cells[1] = (FileCell) FileCellFactory.create(ArtNodeModel.INTERNAL_OUTPUT_PATH + chr + "_F_rev.fq");
    				row = new DefaultRow(key, cells);
    	    	    container.addRowToTable(row);
    				break;
    			case 2:
    	            key = new RowKey("Row " + (row_idx++) + ". " + chr + "_C");
//    	            System.out.println("key 2: "+ key.toString());
    				cells[0] = (FileCell) FileCellFactory.create(ArtNodeModel.INTERNAL_OUTPUT_PATH + chr + "_C_fwd.fq");
    				cells[1] = (FileCell) FileCellFactory.create(ArtNodeModel.INTERNAL_OUTPUT_PATH + chr + "_C_rev.fq");
    				row = new DefaultRow(key, cells);
    	    	    container.addRowToTable(row);
    				break;
    			
    			}
//    	        row_idx++;

    		}
    			

    	    /*DataRow row = new DefaultRow(key, cells);
    	    container.addRowToTable(row);*/
    	    
    	    // check if the execution monitor was canceled
    	    exec.checkCanceled();   

        }
    	// once we are done, we close the container and return its table
        container.close();
    	BufferedDataTable out = container.getTable();
        return new BufferedDataTable[]{out};


	}

    /**
     * Helper method for executing commands
     * @param exec
     * @param map
     * @throws CanceledExecutionException
     * @throws InterruptedException
     * @throws ExecutionException
     * @throws UnsuccessfulExecutionException
     * @throws IOException 
     */
    private void execute_help(final ExecutionContext exec, ArrayList<Art_object> art_obj)
			throws CanceledExecutionException, InterruptedException,
			ExecutionException, UnsuccessfulExecutionException, IOException {
		
		/**
    	 * "$art -i /".$_[$i].".fa -o /".$_[$i]."_ -l 50 -m 200 -f 20 -s 10 -na -p";
    	 */
//    	System.out.println("I CAME TO EXECUTE HELP");
		
		for (int i = 0; i < art_obj.size(); i++) { // #chr ids

			for(int j = 0; j < 6; j++) { // # the trios f0, f1, m0, m1, c0, c1
				ArrayList<String> command = new ArrayList<String>();

				command.add(ArtNodeModel.ART_PATH);
		    	command.add("-na"); // no aln file
		    	command.add("-p"); // enable paired-ended reads

		    	if(m_use_NO_MASK.getBooleanValue()) {
		    		command.add("-nf 0");
		    	}
		    	if(m_use_MATE_PAIR.getBooleanValue()) {
		        	command.add("-mp");
		    	}
		    	if(m_use_ERROR_FREE.getBooleanValue()){
		    		command.add("-ef");
		    	}
		    	if(m_use_SEPERATE_PROFILE.getBooleanValue()){
		    		command.add("-sp");
		    	}
		    	if(m_use_QUIET.getBooleanValue()){
		    		command.add("-q");
		    	}
		    	
		    	command.add("-i");
		    	command.add(art_obj.get(i).getTrio_fasta()[j]); // the name of the fasta file

		    	/**
		    	 * Getting the correct output name: chr21_F_0.fa => chr21_F_0_1/2.fq
		    	 *
		    	 */
		    	String[] dir = art_obj.get(i).getTrio_fasta()[j].split("/");
		    	String prefix = dir[dir.length-1].substring(0, dir[dir.length-1].length()-3); 
		    	
		    	command.add("-o");
		    	command.add(ArtNodeModel.INTERNAL_OUTPUT_PATH + prefix + "_");
		    	
		    	command.add("-l");
		    	command.add(m_LENGTH.getIntValue() + "");
		    	
		    	command.add("-m");
		    	command.add(m_MEAN_SIZE.getIntValue() + "");

		    	command.add("-f");
		    	command.add(m_FOLD.getIntValue() + "");

		    	command.add("-s");
		    	command.add(m_SD.getIntValue() + "");
		    	
		    	command.add("-qs 10");
		    	command.add("-qs2 10");
		    	
//		    	for (String s : command) {
//		    		System.out.print(s + " ");
//		    	}
//		    	System.out.println();
		    	/**Execute**/
//				long startTime = System.currentTimeMillis();
		    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
//				long endTime   = System.currentTimeMillis();
//				DecimalFormat formatter = new DecimalFormat("#0.00000");
//				System.out.println("Execution time for Art " + formatter.format((endTime - startTime) / 1000d) + " seconds");
			}
			mergeFQs(exec, art_obj.get(i).getId());
		}
	}

	
    /**
	 * @param path2id: ids.txt file
	 * @return the ids are retrived from the ids.txt
	 */
	private ArrayList<String> getIds(final String path2id) { 
		
		ArrayList<String> ids = new ArrayList<String>();
		/**
		 * or just ids = FrostRunner.ID_List;
		 */

		try {
			// System.out.println(path2id);
			BufferedReader in = new BufferedReader(new FileReader(path2id));
			String currentLine = "";
			

			while ((currentLine = in.readLine()) != null) {
				if (!currentLine.trim().equals("")) {
					String [] info = currentLine.split("\t");
					int chunks = (Integer.parseInt(info[1])/FrostsLilHelperNodeModel.chunk_length)+1;
					for (int i = 0; i < chunks; i++)
						ids.add(i + "_" + info[0]);
				}
			}
			
		} catch (IOException e) {
			System.out.println("Error when reading " + path2id);
			e.printStackTrace();

		}

		return ids;
	}

	/**
	 * User upload
	 * @param path: The path to the folder where we apparently have the Frost output like it should be:
	 * The fasta files, and the ids.txt
	 * @return
	 */
	private ArrayList<Art_object> getUserInput(final String path, ArrayList<String> ids) {
		// TODO Auto-generated method stub
		ArrayList<Art_object> arrList = new ArrayList<Art_object>();

		/**
		 * Build a Knime datatyp for the trios (TODO)
		 */		
		String directory = path/*.substring(0, path.length()-14)*/;
		File folder = new File(directory); //../../ids.txt
//		System.out.println(folder.toString());
		for(String id : ids) {
//			System.out.println("ID: " + id);
			String[] trios = new String[6];
			for (String fileEntry: folder.list()) {
//				System.out.println(fileEntry.substring(0, fileEntry.length()-8));
				if (fileEntry.endsWith(".fa") && fileEntry.startsWith(id)) {
					trios[0] = directory + "/" + id + "_M0.fa";
					trios[1] = directory + "/" + id + "_M1.fa";
					trios[2] = directory + "/" + id + "_F0.fa";
					trios[3] = directory + "/" + id + "_F1.fa";
					trios[4] = directory + "/" + id + "_C0.fa";
					trios[5] = directory + "/" + id + "_C1.fa";

				}
			}
//			for (String s: trios)
//				System.out.println(s);
			arrList.add(new Art_object(id, trios));
		}
//		for (Art_object a: arrList) {
//			System.out.println(a.getId() + "\t" + a.getTrio_fasta()[0]+ "\t" + a.getTrio_fasta()[1]
//					+ "\t" + a.getTrio_fasta()[2]+ "\t" + a.getTrio_fasta()[3]
//							+ "\t" + a.getTrio_fasta()[4]+ "\t" + a.getTrio_fasta()[5]);
//		}
		
		return arrList;
	}
	
	/**
	 * art_object : retrieve the path to the trio fq file (one id)
	 * @param id 
	 * @throws UnsuccessfulExecutionException 
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 * @throws CanceledExecutionException 
	 * @throws IOException 
	 */
	private void mergeFQs(final ExecutionContext exec, String id) throws CanceledExecutionException, InterruptedException,
	ExecutionException, UnsuccessfulExecutionException, IOException {
			
		for (int i = 0; i < 3; i++) {
			String indiv = "";
			switch(i) {
			case 0:
				indiv = "M";
				break;
			case 1:
				indiv = "F";
				break;			
			case 2:
				indiv = "C";
				break;
			}
			for (int j = 0; j < 2; j++) {
				String read = ""; int num = 0;
				switch(j) {
				case 0:
					read = "fwd";
					num = 1;
					break;
				case 1:
					read = "rev";
					num = 2;
					break;
				}
				String[] input_fqs = {ArtNodeModel.INTERNAL_OUTPUT_PATH + id + "_" + indiv + "0_" + num + ".fq",
						ArtNodeModel.INTERNAL_OUTPUT_PATH + id + "_" + indiv + "1_" + num + ".fq"};
				String output_fq = ArtNodeModel.INTERNAL_OUTPUT_PATH + id + "_" + indiv + "_" + read + ".fq";
				
				ArrayList<String> merge_command = new ArrayList<String>(5);
				merge_command.add("sh /home/ibis/tanzeem.haque/Documents/Scripts/Sequenciator/mergeFqs.sh");
				merge_command.add(input_fqs[0]);
				merge_command.add(input_fqs[1]);
				merge_command.add(output_fq);
				
//				for (String s : merge_command) {
//					System.out.print(s + " ");
//				}
//				System.out.println();
				
				
				/**Execute**/
		    	Executor.executeCommand(new String[]{StringUtils.join(merge_command, " ")},exec,LOGGER);
		    	boolean fq1_del = new File(input_fqs[0]).delete(), fq2_del = new File(input_fqs[1]).delete();
		    	System.out.println("FQDEL: " + fq1_del + "\t" + fq2_del);
		    	
			}
		}
	}
	
    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO Code executed on reset.
        // Models build during execute are cleared here.
        // Also data handled in load/saveInternals will be erased here.
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
        
        // TODO: check if user settings are available, fit to the incoming
        // table structure, and the incoming types are feasible for the node
        // to execute. If the node can execute in its current state return
        // the spec of its output data table(s) (if you can, otherwise an array
        // with null elements), or throw an exception with a useful user message
    	// flowVariables
        
    	
    	/**
    	 * Check OptionalInputPort
    	 */
		try{
			inSpecs[0].getColumnNames();
			String[] colNames = inSpecs[0].getColumnNames();
			optionalPort=true;
			
//			for (String s : colNames)

			if (colNames.length==6 && colNames[0].equals(FrostNodeModel.OUT_COL + "M0.fa") &&
					colNames[5].equals(FrostNodeModel.OUT_COL + "C1.fa")) 
				m_ID_PATH.setEnabled(false);
			else
				throw new InvalidSettingsException("The in-port is invalid! Don't use an in-port or connect the right node.");
			

    	}catch(NullPointerException npe){
    		m_ID_PATH.setEnabled(true);
			optionalPort=false;
    	
    	}
		System.out.println("Optional Port: " + optionalPort + "\t" + "Enabled: " + m_ID_PATH.isEnabled());

		int col_num = 2; //ArtNodeModel.IDS.size();
    	DataColumnSpec[] allColSpecs = new DataColumnSpec[col_num];/**
    	 * create the column(s) for (multi)fastq output of a trio with 6 rows
    	 */
    	/**
    	for (int i = 0; i < col_num; i++) {
    	    allColSpecs[i] = new DataColumnSpecCreator("Col. " + (i+1) + ": " + OUT_COL + ArtNodeModel.IDS.get(i), FileCell.TYPE).createSpec();

    	}**/
    	allColSpecs[0] = new DataColumnSpecCreator("Col. 1: " + OUT_COL + "fwd", FileCell.TYPE).createSpec();
    	allColSpecs[1] = new DataColumnSpecCreator("Col. 2: " + OUT_COL + "rev", FileCell.TYPE).createSpec();

    	DataTableSpec outputSpec = new DataTableSpec(allColSpecs);
//		System.out.println(optionalPort);
        return new DataTableSpec[]{outputSpec};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {

        // TODO save user settings to the config object.
        if (!optionalPort)
        	m_ID_PATH.saveSettingsTo(settings);
        m_LENGTH.saveSettingsTo(settings);
        m_MEAN_SIZE.saveSettingsTo(settings);
        m_FOLD.saveSettingsTo(settings);
        m_SD.saveSettingsTo(settings);
        m_use_NO_MASK.saveSettingsTo(settings);
        m_use_MATE_PAIR.saveSettingsTo(settings);
        m_use_ERROR_FREE.saveSettingsTo(settings);
        m_use_SEPERATE_PROFILE.saveSettingsTo(settings);
        m_use_QUIET.saveSettingsTo(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
        // TODO load (valid) settings from the config object.
        // It can be safely assumed that the settings are valided by the 
        // method below.
        
    	if (!optionalPort)
        	m_ID_PATH.loadSettingsFrom(settings);
        m_LENGTH.loadSettingsFrom(settings);
        m_MEAN_SIZE.loadSettingsFrom(settings);
        m_FOLD.loadSettingsFrom(settings);
        m_SD.loadSettingsFrom(settings);
        m_use_NO_MASK.loadSettingsFrom(settings);
        m_use_MATE_PAIR.loadSettingsFrom(settings);
        m_use_ERROR_FREE.loadSettingsFrom(settings);
        m_use_SEPERATE_PROFILE.loadSettingsFrom(settings);
        m_use_QUIET.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
        // TODO check if the settings could be applied to our model
        // e.g. if the count is in a certain range (which is ensured by the
        // SettingsModel).
        // Do not actually set any values of any member variables.
    	if (!optionalPort)
        	m_ID_PATH.validateSettings(settings);
        m_LENGTH.validateSettings(settings);
        m_MEAN_SIZE.validateSettings(settings);
        m_FOLD.validateSettings(settings);
        m_SD.validateSettings(settings);
        m_use_NO_MASK.validateSettings(settings);
        m_use_MATE_PAIR.validateSettings(settings);
        m_use_ERROR_FREE.validateSettings(settings);
        m_use_SEPERATE_PROFILE.validateSettings(settings);
        m_use_QUIET.validateSettings(settings);

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        
        // TODO load internal data. 
        // Everything handed to output ports is loaded automatically (data
        // returned by the execute method, models loaded in loadModelContent,
        // and user settings set through loadSettingsFrom - is all taken care 
        // of). Load here only the other internals that need to be restored
        // (e.g. data used by the views).

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
       
        // TODO save internal models. 
        // Everything written to output ports is saved automatically (data
        // returned by the execute method, models saved in the saveModelContent,
        // and user settings saved through saveSettingsTo - is all taken care 
        // of). Save here only the other internals that need to be preserved
        // (e.g. data used by the views).

    }
    
    /**
     * 
     * A nested internal class for the art inputs
     * the chromosome id along with the 6 trio fastas
     * @author tanzeem.haque
     *
     */
    private static class Art_object {
    	private String id;
    	private String[] trio_fasta = new String[6];
    	
    	private Art_object(String id, String[] trio_fasta) {
    		setId(id);
    		setTrio_fasta(trio_fasta);
    	}

		public String getId() {
			return id;
		}

		public void setId(String id) {
			this.id = id;
		}

		public String[] getTrio_fasta() {
			return trio_fasta;
		}

		public void setTrio_fasta(String[] trio_fasta) {
			this.trio_fasta = trio_fasta;
		}
    }

}

package de.helmholtz_muenchen.ibis.ngs.fastSam2Bam;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.SettingsStorageNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.samsplitter.helpers.FileHelpers;
import de.helmholtz_muenchen.ibis.utils.ngs.samsplitter.samhelpers.SamSplitter;
import de.helmholtz_muenchen.ibis.utils.threads.ExecuteThread;
import picard.sam.MergeSamFiles;

/**
 * This node is a faster implementation of a sam to bam converter based on
 * a script from Jonathan Hoser. 
 *
 * @author Michael Kluge
 */
public class FastSam2BamNodeModel extends SettingsStorageNodeModel {
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FastSam2BamNodeModel.class);
	
	// keys for the settings models
	protected static final String CFGKEY_GENOME_FILE		= "FS2B_genomeFile";
	protected static final String CFGKEY_OUTPUT_PATH 		= "FS2B_outputPath";
	protected static final String CFGKEY_PATH_SAMTOOLS 		= "FS2B_samTools";
//	protected static final String CFGKEY_PATH_PICTOOLS 		= "FS2B_picTools";
	protected static final String CFGKEY_TMP_PATH 			= "FS2B_tmpPath";
	protected static final String CFGKEY_USE_RAM 			= "FS2B_useRam";
	protected static final String CFGKEY_DELETE_SAM 		= "FS2B_deleteSam";
	protected static final String CFGKEY_CORE_NUMBER 		= "FS2B_threadNumber";
	protected static final String CFGKEY_SPLIT_SIZE 		= "FS2B_splitSize";
	
	// default values
	protected static final String DEFAULT_GENOME_FILE 		= "";
	protected static final String DEFAULT_OUTPUT_PATH 		= "";
	protected static final String DEFAULT_PATH_SAMTOOLS 	= "";
	protected static final String DEFAULT_PATH_PICTOOLS 	= "";
	protected static final String DEFAULT_TMP_PATH 			= "/tmp";
	protected static final boolean DEFAULT_USE_RAM 			= false;
	protected static final boolean DEFAULT_DELETE_SAM 		= false;
	protected static final int DEFAULT_CORE_NUMBER 			= 2;
	protected static final int DEFAULT_SPLIT_SIZE 			= 1000000;
	protected static final String DEFAULT_USE_RAM_PATH 		= "/dev/shm/";
	
	// other settings
//	private final static String SAMTOOLS_BINARY_NAME = "samtools";
//	private final static String PICTOOLS_JAR_NAME = "MergeSamFiles.jar";
	
	// name of the output variables
	public static final String OUT_COL1 = "Path2BAMFile";
	public static final String OUT_COL2 = "Path2Bai";
	public static final String OUT_COL3 = "Path2SEQFile";
       
	// definition of SettingsModel (all prefixed with SET)
	private final SettingsModelString SET_GENOME 			= new SettingsModelString(CFGKEY_GENOME_FILE, DEFAULT_GENOME_FILE);
    private final SettingsModelString SET_OUTPUT_PATH		= new SettingsModelString(CFGKEY_OUTPUT_PATH, DEFAULT_OUTPUT_PATH);
    private final SettingsModelInteger SET_CORE_NUMBER		= new SettingsModelInteger(CFGKEY_CORE_NUMBER, DEFAULT_CORE_NUMBER);
    private final SettingsModelInteger SET_SPLIT_SIZE		= new SettingsModelInteger(CFGKEY_SPLIT_SIZE, DEFAULT_SPLIT_SIZE);
    private final SettingsModelString SET_TMP_PATH			= new SettingsModelString(CFGKEY_TMP_PATH, DEFAULT_TMP_PATH);
    private final SettingsModelBoolean SET_USE_RAM_AS_TMP  	= new SettingsModelBoolean(CFGKEY_USE_RAM, DEFAULT_USE_RAM);
    private final SettingsModelString SET_PATH_SAMTOOLS		= new SettingsModelString(CFGKEY_PATH_SAMTOOLS, DEFAULT_PATH_SAMTOOLS);
//    private final SettingsModelString SET_PATH_PICTOOLS    	= new SettingsModelString(CFGKEY_PATH_PICTOOLS, DEFAULT_PATH_PICTOOLS);
    private final SettingsModelBoolean SET_DELETE_SAM		= new SettingsModelBoolean(CFGKEY_DELETE_SAM, DEFAULT_DELETE_SAM);

    
    /**
     * Constructor for the node model.
     */
    protected FastSam2BamNodeModel() {
        super(1, 1);
        addSetting(SET_GENOME);
    	addSetting(SET_OUTPUT_PATH);
    	addSetting(SET_CORE_NUMBER);
    	addSetting(SET_SPLIT_SIZE);
    	addSetting(SET_TMP_PATH);
    	addSetting(SET_USE_RAM_AS_TMP);
    	addSetting(SET_PATH_SAMTOOLS);
//    	addSetting(SET_PATH_PICTOOLS);
    	addSetting(SET_DELETE_SAM);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception {
    	exec.setProgress(0.01); // tell the user that we started with the work
    	int cores = SET_CORE_NUMBER.getIntValue();
    	String nameOfBamFile = "";
    	String nameOfBaiFile = "";
    	
    	// create tmp folder
    	String randomID = this.getClass().getSimpleName() + "_" + (System.currentTimeMillis() / 1000L) + "_" + new Random().nextInt(10000);
    	File tmpPath = new File(SET_TMP_PATH.getStringValue() + File.separator + randomID);
    	if(!tmpPath.mkdirs()) 
    		throw new FileNotFoundException("Could not create temp path '" + tmpPath.getAbsolutePath() + "'.");
    	
    	// get input file
    	Iterator<DataRow> it = inData[0].iterator();
    	if(it.hasNext()) {
			// create worker pool
			ExecutorService pool = Executors.newFixedThreadPool(cores);			
    		String inputFile = it.next().getCell(0).toString();
    		String outfolder = SET_OUTPUT_PATH.getStringValue();
    		if(outfolder.equals("") || Files.notExists(Paths.get(outfolder))) {
    			outfolder = new File(inputFile).getParent();
    		}
    		
        	nameOfBamFile = outfolder + File.separator + new File(inputFile).getName() + ".sorted.bam";
        	nameOfBaiFile = nameOfBamFile.replaceFirst(".bam$", ".bai");
        	
        	
    		// check if input file is there
    		if(!new File(inputFile).exists()) 
    			throw new FileNotFoundException("Input BAM File '" + inputFile + "' not found.");
    		
    		// check, if bam file is already there
    		if(!(new File(nameOfBamFile).exists() && new File(nameOfBaiFile).exists())) {
    			new File(outfolder).mkdirs();
	    		// run SamSplitter
    			FileHelpers.noExitRule = true; 	// do not call System.exit()!
    			SamSplitter.noExit = true; 		// do not call System.exit()!
				SamSplitter ssr = new SamSplitter();
				ssr.setInFile(inputFile);
				ssr.setOutDir(tmpPath.getAbsolutePath());
				ssr.setSplitSize(SET_SPLIT_SIZE.getIntValue());
				ssr.split();
				exec.setProgress(0.1);
				
				ArrayList<String> partsBam = new ArrayList<String>();
				ArrayList<String> partsBamSorted = new ArrayList<String>();
				ArrayList<ExecuteThread> tasks = new ArrayList<ExecuteThread>();
				
				////////////////////////// Convert SAM -> BAM ///////////////////////////////////////////////
				// get file of folder which must be processed
				ArrayList<String> partsSam = IO.getFilesOfFolder(tmpPath.getAbsolutePath());			
				for(String parFile : partsSam) {
//					String command = SET_PATH_SAMTOOLS.getStringValue() + File.separator;
					System.out.println(parFile);
					
					String bamFile = parFile + ".bam";
					partsBam.add(bamFile);
					String command = SET_PATH_SAMTOOLS.getStringValue()+ " view -b -S -T " + SET_GENOME.getStringValue() + " -o " + bamFile + " " + parFile;
					String outFile = parFile + ".log";
					ExecuteThread task = new ExecuteThread(new String[] {command}, LOGGER, outFile, outFile, null, null, null, null,null);
					pool.submit(task);
					tasks.add(task);
				}
				pool.shutdown();
				// wait until all tasks were executed
				while (!pool.awaitTermination(5, TimeUnit.SECONDS)) { 
					 // check, if operation was canceled			 
					 try {
						 exec.checkCanceled();
					 } catch (CanceledExecutionException e){
						// kill jobs
						for(ExecuteThread task : tasks)
							task.cancel();
						
						pool.shutdownNow();
						while (!pool.isTerminated()) {
						}
						throw(new CanceledExecutionException("Execution of FastSam2Bam node was canceled."));
					}
					Thread.sleep(1000);
				}
				exec.setProgress(0.2);
				/////////////////////////////////////////////////////////////////////////////////////////////
				
				// delete the sam split files
				tasks.clear();
				for(String parFile : partsSam) {
					new File(parFile).delete();
				}
				
				//////////////////////////// Sort BAM ////////////////////////////////////////////////////////
				// re-init worker pool
				pool = Executors.newFixedThreadPool(cores);
				for(String bamFile : partsBam) {
//					String command = SET_PATH_SAMTOOLS.getStringValue() + File.separator;
					String sortedBamFile = bamFile + ".sorted.bam";
					partsBamSorted.add(sortedBamFile);
					String command = SET_PATH_SAMTOOLS.getStringValue() + " sort -m 10000000000 -o" + sortedBamFile + " " + bamFile;
					String outFile = partsBam + ".log";
					ExecuteThread task = new ExecuteThread(new String[] {command}, LOGGER, outFile, outFile, null, null, null, null,null);
					pool.submit(task);
					tasks.add(task);
				}
				pool.shutdown();
				// wait until all tasks were executed
				while (!pool.awaitTermination(5, TimeUnit.SECONDS)) { 
					 // check, if operation was canceled			 
					 try {
						 exec.checkCanceled();
					 } catch (CanceledExecutionException e){
						// kill jobs
						for(ExecuteThread task : tasks)
							task.cancel();
						
						pool.shutdownNow();
						while (!pool.isTerminated()) {
						}
						throw(new CanceledExecutionException("Execution of FastSam2Bam node was canceled."));
					}
					Thread.sleep(1000);
				}
				exec.setProgress(0.3);
		    	/////////////////////////////////////////////////////////////////////////////////////////////
		    	
		    	// delete the bam split files
				tasks.clear();
				for(String parFile : partsBam) {
					new File(parFile).delete();
				}
				
				////////////////////////////Merge it again ////////////////////////////////////////////////////////
				int parts = partsBamSorted.size();
				String [] args = new String[7+parts];
				for(int i = 0; i< parts; i++) {
					args[i] = "INPUT=" + partsBamSorted.get(i);
				}
				
				args[parts+0] = "OUTPUT=" + nameOfBamFile;
				args[parts+1] = "SORT_ORDER=coordinate";
				args[parts+2] = "USE_THREADING=true";
				args[parts+3] = "CREATE_INDEX=true";
				args[parts+4] = "MAX_RECORDS_IN_RAM=600000000";
				args[parts+5] = "TMP_DIR=" + tmpPath;
				args[parts+6] = "VALIDATION_STRINGENCY=LENIENT";
				new MergeSamFiles().instanceMain(args);
				
//				ArrayList<String> mergeCommand = new ArrayList<String>();
//				mergeCommand.add("java -jar");
//				mergeCommand.add(SET_PATH_PICTOOLS.getStringValue() + File.separator + PICTOOLS_JAR_NAME);
//				for(String bamSorted : partsBamSorted) {
//					mergeCommand.add("INPUT=" + bamSorted + ".bam");
//				}
//				mergeCommand.add("OUTPUT=" + nameOfBamFile);
//				mergeCommand.add("SORT_ORDER=coordinate");
//				mergeCommand.add("USE_THREADING=true");
//				mergeCommand.add("CREATE_INDEX=true");
//				mergeCommand.add("MAX_RECORDS_IN_RAM=600000000");
//				mergeCommand.add("TMP_DIR=" + tmpPath);
//				mergeCommand.add("VALIDATION_STRINGENCY=LENIENT");
//				String outLogFile = nameOfBamFile + ".log";
//				Executor.executeCommand(new String[]{StringUtils.join(mergeCommand, " ")}, exec, null, LOGGER, outLogFile, outLogFile, null, null);
    		}
    		
	    	// ensure that bam file is there
			File bamFile = new File(nameOfBamFile);
			if(bamFile.exists() && bamFile.getTotalSpace() > 1024) { // file is not empty
		    	// check, if input file should be deleted
		    	if(this.SET_DELETE_SAM.getBooleanValue())
		    		new File(inputFile).delete();
			}
			else {
				throw new IllegalArgumentException("Bam file could not be sucessfully created. Please check the log files in the temp folder...");
			}
	    }
    	
    	
    	exec.setProgress(0.95);
    	// delete the temp stuff
    	ArrayList<String> tmpFiles = IO.getFilesOfFolder(tmpPath.getAbsolutePath());	
    	for(String file : tmpFiles) 
    		new File(file).delete();
    	tmpPath.delete();
    	///////////////////////////////////////////////////////
    	
    	exec.setProgress(1); // we are done :)

    	// write output
    	BufferedDataContainer cont= exec.createDataContainer(getDataOutSpec1());
    	FileCell[] c = new FileCell[]{ (FileCell)FileCellFactory.create(nameOfBamFile), (FileCell)FileCellFactory.create(nameOfBaiFile), (FileCell)FileCellFactory.create(SET_GENOME.getStringValue()) };
    	cont.addRowToTable(new DefaultRow("Row0", c));
    	cont.close();
        return new BufferedDataTable[]{cont.getTable()};
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
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	// validate binaries
//    	IO.isBinaryValid(this.SET_PATH_SAMTOOLS.getStringValue() + File.separator + SAMTOOLS_BINARY_NAME, true);
    	IO.isBinaryValid(this.SET_PATH_SAMTOOLS.getStringValue(), true);

    	//    	IO.isBinaryValid(this.SET_PATH_PICTOOLS.getStringValue() + File.separator + PICTOOLS_JAR_NAME, false);

    	// TODO: check fasta genome file 
    	// TODO: check other stuff ?!?
    	
        return new DataTableSpec[]{getDataOutSpec1()};
    }
    
    /**
     * returns the first output specifications of this node
     * @return
     */
    private DataTableSpec getDataOutSpec1() {
    	return new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, BAMCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL3, FileCell.TYPE).createSpec()});
    }

	@Override
	protected void loadInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {
	}

	@Override
	protected void saveInternals(File nodeInternDir, ExecutionMonitor exec) throws IOException, CanceledExecutionException {	
	}
}


package de.helmholtz_muenchen.ibis.ngs.bwa;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

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
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.BinaryHandler;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;



/**
 * This is the model implementation of BWA.
 * 
 * @author Jan Quell
 * @author Maximilian Hastreiter
 */
public class BWANodeModel extends HTExecutorNodeModel {
    
	public static final String CFGKEY_USEPREFPAGE = "usePrefPage";
	public static final String CFGKEY_REFSEQFILE = "refseqfile";
	public static final String CFGKEY_BWAFILE = "bwafile";
	public static final String CFGKEY_CHECKCOLORSPACED = "checkColorSpaced";
	public static final String CFGKEY_BWTINDEX = "bwtIndex";
	public static final String CFGKEY_READTYPE = "readType";
	public static final String CFGKEY_CHECKINDEX = "checkIndexRefSeq";
	public static final String CFGKEY_READGROUP = "readgroup";
	public static final String CFGKEY_READGROUPBOOLEAN = "readgroupboolean";
	public static final String CFGKEY_ALNALGO = "alnalgo";
	public static final String CFGKEY_THREADS = "alnthreads";
	
	
    // definition of SettingsModel
	private final SettingsModelString m_refseqfile = new SettingsModelString(CFGKEY_REFSEQFILE,"");
	private final SettingsModelString m_bwafile = new SettingsModelString(CFGKEY_BWAFILE,"");
	private final SettingsModelBoolean m_checkIndexRefSeq = new SettingsModelBoolean(CFGKEY_CHECKINDEX,true);
	private final SettingsModelBoolean m_checkColorSpaced = new SettingsModelBoolean(CFGKEY_CHECKCOLORSPACED, false);
	private final SettingsModelString m_bwtIndex = new SettingsModelString(CFGKEY_BWTINDEX,"BWT-SW");
	private final SettingsModelString m_alnalgo = new SettingsModelString(CFGKEY_ALNALGO,"BWA-MEM");
	private final SettingsModelString m_readType = new SettingsModelString(CFGKEY_READTYPE,"auto-detect");
	private final SettingsModelString m_readGroup = new SettingsModelString(CFGKEY_READGROUP,"@RG\\tID:foo\\tSM:bar");
	private final SettingsModelBoolean m_readGroupBoolean = new SettingsModelBoolean(CFGKEY_READGROUPBOOLEAN,false);
	private final SettingsModelIntegerBounded m_ALN_THREADS = new SettingsModelIntegerBounded(CFGKEY_THREADS,4, 1, Integer.MAX_VALUE);
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(BWANodeModel.class);
	
	//The Output Col Names
	public static final String OUT_COL1 = "Path2SAMFile";
	public static final String OUT_COL2 = "Path2RefFile";
	

	
	
    /**
     * Constructor for the node model.
     */
    protected BWANodeModel() {
    	
    	super(1, 1);
    	
    	m_readType.setEnabled(false);
    	
    }

    static SettingsModelString createSettingsModelSelection() {
    	return new SettingsModelString("bwa-path","");
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	    
    	/**Initialize logfile**/
    	String fle = inData[0].iterator().next().getCell(0).toString();
    	String logfile = fle.substring(0,fle.lastIndexOf("/")+1)+"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("BWA"));
    	/**end initializing logfile**/
    	
    	/**
		 * Get the Parameters
		 */
    	String path2refFile = m_refseqfile.getStringValue();
    	String path2readFile = inData[0].iterator().next().getCell(0).toString();
    	String path2readFile2 = inData[0].iterator().next().getCell(1).toString();
    	String readTypePrevious = getAvailableInputFlowVariables().get("readType").getStringValue();
    	String readType = m_readType.getStringValue();
    	
    	String basePath = path2readFile.substring(0,path2readFile.lastIndexOf('/')+1);
    	String outBaseName1 = path2readFile.substring(path2readFile.lastIndexOf("/")+1,path2readFile.lastIndexOf("."));
    	String outBaseName = outBaseName1;
    	String outBaseName2 = outBaseName1;
    	String memOut = basePath+outBaseName1+"_mem.sam";
    	
    	if(path2readFile2.length() > 1 && !path2readFile2.equals("na")) {
    		outBaseName2 = path2readFile2.substring(path2readFile2.lastIndexOf("/")+1,path2readFile2.lastIndexOf("."));
	    	if(!path2readFile.equals(path2readFile2)) {
	    		outBaseName = outBaseName1 + "_" + outBaseName2;
	    	}
    	}else{
    		
    	}
    	
    	String out2Name = basePath+outBaseName+"_aln.sam";
    	String out1Name = basePath+outBaseName+"_aln_sa.sai";
    	String out11Name = basePath+outBaseName1+"_aln_sa_1.sai";
    	String out12Name = basePath+outBaseName2+"_aln_sa_2.sai";

    	Boolean isBam = false;
    	    	
    	if(path2readFile.substring(path2readFile.length()-3, path2readFile.length()) == "bam") {isBam = true;}
    	if(!readTypePrevious.equals("") && !readTypePrevious.equals(readType)) {readType = readTypePrevious;}
    	if(isBam) {path2readFile2 = path2readFile;}
    	
    	String path2bwa = m_bwafile.getStringValue();
    	int threads = m_ALN_THREADS.getIntValue();
    	
    	if(path2readFile2.length() > 1 && !path2readFile2.equals("na")) {
    		outBaseName2 = path2readFile2.substring(path2readFile2.lastIndexOf("/")+1,path2readFile2.lastIndexOf("."));
	    	if(!path2readFile.equals(path2readFile2)) {
	    		outBaseName = outBaseName1 + "_" + outBaseName2;
	    	}
    	}
    	String colorSpaced = "-c ";

    	
    	if(path2readFile.substring(path2readFile.length()-3, path2readFile.length()) == "bam") {
    		isBam = true;
    	}
    	
    	if(!readTypePrevious.equals("") && !readTypePrevious.equals(readType)) {
    		readType = readTypePrevious;
    	}
    	LOGGER.info("Read Type: " + readType + "\n");
    	/*******************************************************************************************/   	
    	
    	/**
    	 * Create Index for Read Files
    	 */

    	bwa_index(exec,logBuffer, colorSpaced, path2bwa, path2refFile, path2readFile2);

    	/**
    	 * Run bwa aln
    	 */
    	
    	
    	if(!m_alnalgo.getStringValue().equals("BWA-MEM")){
        	LOGGER.info("Find the SA coordinates of the input reads.\n");
        	bwa_aln(exec,readType, basePath, outBaseName, outBaseName1, outBaseName2, path2refFile, path2bwa, path2readFile, logBuffer, path2readFile2, isBam,threads);
        	LOGGER.info("Finished BWA aln...");
    	}

    	
    	
    	/**
    	 * Map Reads
    	 */
    	bwa_map(exec,readType,logBuffer,path2bwa,path2refFile,path2readFile,out1Name,out2Name,out11Name,out12Name,path2readFile2,memOut);
    	//TODO CHECK bwa_map method for errors !!!! TODO
    	
    	/**
    	 * OUTPUT
    	 */
    	if(m_alnalgo.getStringValue().equals("BWA-MEM")){
    		out2Name = memOut;
    	}
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(out2Name),
    			(FileCell) FileCellFactory.create(path2refFile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	pushFlowVariableString("BAMSAMINFILE",out2Name);

		return new BufferedDataTable[]{outTable};
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
    }

    /**
     * Runs bwa index
     * @param logBuffer
     * @param colorSpaced
     * @param path2bwa
     * @param path2refFile
     * @param path2readFile
     * @throws Exception 
     */
    private void bwa_index(ExecutionContext exec,StringBuffer logBuffer, String colorSpaced, String path2bwa, String path2refFile, String path2readFile) throws Exception{
    	/**Only execute if Index needs to be created**/
    	if(m_checkIndexRefSeq.getBooleanValue()) {
    		
    		LOGGER.info("Indexing reference sequence.\n");
    		
        	ArrayList<String> command = new ArrayList<String>();
        	// Constant values
        	command.add(path2bwa+" index");
        	
        	//Indexing Type
	    	if(m_bwtIndex.getStringValue() == "BWT-SW") {
	    		command.add("-a bwtsw");
	    	} else {
	    		command.add("-a is");
	    	}
	    	// Colorspace
	    	if(m_checkColorSpaced.getBooleanValue()) {
	    		command.add(colorSpaced);
	    	}
	    	
	    	//Add Reference genome
	    	command.add(path2refFile);

	    	/**Execute**/
	    	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, null, null, null, null, null, null, null);
//	    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
			
    	} else {
    		LOGGER.info("Indexing reference sequence SKIPPED.\n");
    	}
    }
    
    private void bwa_aln(ExecutionContext exec, String readType, String basePath, String outBaseName, String outBaseName1, String outBaseName2, String path2refFile, String path2bwa, String path2readFile, StringBuffer logBuffer, String path2readFile2, boolean isBam, int Threads) throws Exception{
    	
    	
    	ArrayList<String> command = new ArrayList<String>();
    	// Constant values
    	command.add(path2bwa+" aln");
    	
    	String outName = basePath+outBaseName+"_aln_sa.sai";
    	String out11Name = basePath+outBaseName1+"_aln_sa_1.sai";
    	String out12Name = basePath+outBaseName2+"_aln_sa_2.sai";
    	String outfile = outName;
    	
    	//Multi-Threading 
    	command.add("-t " +Threads);
    	
    	//If Inputfile is in bam format
    	if(readType.equals("paired-end")){
    		outfile = out11Name;				//Set Outfile for forward reads
    		if(isBam){
    			command.add("-b1");
    		}
    	}else{	// Single-end
    		if(isBam){
    			command.add("-b0");
    		}
		}
    	

    	
    	//Perform aln for forward reads OR single end reads
    	command.add(path2refFile);
    	command.add(path2readFile);
    	command.add("-f "+outfile);
    	/**Execute**/
    	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, null, null, null, null, null, null, null);
//    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
   
		//If paired end, repeat previous step
    	if(readType.equals("paired-end")) {
    		
        	if(isBam){
        		command.set(2, "-b2");        	
        		command.set(4, path2readFile2);
            	command.set(5, " -f "+ out12Name);
        	}else{
        		command.set(3, path2readFile2);
            	command.set(4, " -f "+ out12Name);
        	}
        	/**Execute**/
        	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, null, null, null, null, null, null, null);
//        	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
		}
    }
    
    
  private void bwa_map(ExecutionContext exec,String readType, StringBuffer logBuffer, String path2bwa, String path2refFile, String path2readFile, String out1Name, String out2Name, String out11Name, String out12Name, String path2readFile2, String memOut) throws Exception{ 	
    	
  		ArrayList<String> command = new ArrayList<String>();
  		String alnalgo = m_alnalgo.getStringValue(); 	
    	
// ### BWA-Backtrack ###
    	if(alnalgo.equals("BWA-backtrack")) {
    		
        	if(readType.equals("single-end")) {
        		// bwa samse sequence.fasta aln_sa.sai s_1_1_sequence.txt > aln.sam
        		LOGGER.info("Generate alignments in the SAM format given single-end reads.\n");

        		command.add(path2bwa+" samse");
        		
        		/**Other Options**/
        		if(m_readGroupBoolean.getBooleanValue()){
        			command.add("-r "+m_readGroup.getStringValue());
        		}
        		/**In and Outfiles**/
        		command.add("-f "+out2Name);
        		command.add(path2refFile);
        		command.add(out1Name);
        		command.add(path2readFile);

        	} else {
        		// bwa sampe sequence.fasta aln_sa_1.sai aln_sa_2.sai s_1_1_sequence.fq s_1_2_sequence.fq > aln.sam
        		LOGGER.info("Generate alignments in the SAM format given paired-end reads.\n");
        		
        		command.add(path2bwa+" sampe");
        		
        		/**Other Options**/
        		if(m_readGroupBoolean.getBooleanValue()){
        			command.add("-r "+m_readGroup.getStringValue());
        		}
        		/**In and Outfiles**/
        		command.add("-f "+out2Name);
        		command.add(path2refFile);
        		command.add(out11Name);
        		command.add(out12Name);
        		command.add(path2readFile);
        		command.add(path2readFile2);
        		
        	}
// ### BWA-SW ###
    	} else if(alnalgo.equals("BWA-SW")) {
        	LOGGER.info("Generate alignments in the SAM format.\n");

        	command.add(path2bwa+" bwasw");
        	
    		/**In and Outfiles**/
        	command.add("-f "+out2Name);
        	command.add(path2refFile);
        	command.add(path2readFile);
    		if(readType.equals("paired-end")) {
    			command.add(path2readFile2);
    		} 		
// ### BWA-MEM ###
    	} else if(alnalgo.equals("BWA-MEM")) {
        	LOGGER.info("Generate alignments in the SAM format.\n");

        	command.add(path2bwa+" mem");
        	command.add("-t "+m_ALN_THREADS.getIntValue());
    		if(m_readGroupBoolean.getBooleanValue()){
    			command.add("-R "+m_readGroup.getStringValue());
    		}
        	
    		/**In and Outfiles**/
        	command.add(path2refFile);
        	command.add(path2readFile);
    		if(readType.equals("paired-end")) {
    			command.add(path2readFile2);
    		} 
    	}
		
    	
    	/** check if run was already sucessful **/
//    	String[] com = command.toArray(new String[command.size()]);
    	File lockFile = new File(path2readFile.substring(0,path2readFile.lastIndexOf(".")) + ".BWA" +  SuccessfulRunChecker.LOCK_ENDING);
//    	String lockCommand = ExecuteThread.getCommand(com);
//    	boolean terminationState = SuccessfulRunChecker.hasTerminatedSuccessfully(lockFile, lockCommand);
//		LOGGER.info("Successful termination state: " + terminationState);
    	
		// do not execute if termination state is true
//		if(!terminationState) {
//			SuccessfulRunChecker checker = new SuccessfulRunChecker(lockFile, lockCommand);
		
	    	/**Execute**/
	    	if(alnalgo.equals("BWA-MEM")) {
//	    		Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER,memOut);
	        	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, null, lockFile, memOut, null, null, null, null);

	    	}else{
//	    		Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
	        	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, null, lockFile, null, null, null, null, null);

	    	}
//	    	checker.writeOK();
		}
//	}
    
    
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	String toolPath = BinaryHandler.checkToolAvailability("bwa");
    	if(toolPath == null) {
    		toolPath = "";
    	}
    	m_bwafile.setStringValue(toolPath);
    	
    	// Warning if there is a problem with readType
    	String readTypePrevious = getAvailableInputFlowVariables().get("readType").getStringValue();
    	String readType = m_readType.getStringValue();
    	if(!readTypePrevious.equals("") && !readTypePrevious.equals(readType) && !readType.equals("auto-detect")) {
    		setWarningMessage("The previous node indicates that you have " + readTypePrevious + " reads, but you have chosen " + readType + ". BWA will use " + readTypePrevious + " mapping.");
    	}
    	
        //Version control
        if(FileValidator.versionControl(m_bwafile.getStringValue(),"BWA")==1){
        	setWarningMessage("WARNING: You are using a newer BWA version than "+FileValidator.BWA_VERSION +"! This may cause problems");
        }else if(FileValidator.versionControl(m_bwafile.getStringValue(),"BWA")==2){
        	throw new InvalidSettingsException("You are using a outdated version of BWA! Please update your version");
        }else if(FileValidator.versionControl(m_bwafile.getStringValue(),"BWA")==-1){
        	System.out.println("Something wrong here");
        }
    	
    	if(m_refseqfile.getStringValue().length() > 1) {
	    	if(!FileValidator.checkFastaFormat(m_refseqfile.getStringValue())){
	            throw new InvalidSettingsException("Reference (genome) sequence file is not in FastA format or does not contain nucleotide sequences!");
	    	}
    	}
    	
    	// Check input ports
    	String[] cn=inSpecs[0].getColumnNames();
    	if(!cn[0].equals("") && !cn[0].equals("Path2ReadFile1")) {
    		throw new InvalidSettingsException("This node is incompatible with the previous node. The outport of the previous node has to fit to the inport of this node.");
    	}

        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	/** added for HTE **/
    	super.saveSettingsTo(settings);
    	
    	m_bwafile.saveSettingsTo(settings);
    	m_refseqfile.saveSettingsTo(settings);
    	m_bwtIndex.saveSettingsTo(settings);
    	m_readType.saveSettingsTo(settings);
    	m_checkColorSpaced.saveSettingsTo(settings);
    	m_checkIndexRefSeq.saveSettingsTo(settings);
    	m_readGroup.saveSettingsTo(settings);
    	m_readGroupBoolean.saveSettingsTo(settings);
    	m_alnalgo.saveSettingsTo(settings);
    	m_ALN_THREADS.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	/** added for HTE **/
    	super.loadValidatedSettingsFrom(settings);
    	
    	m_bwafile.loadSettingsFrom(settings);
    	m_refseqfile.loadSettingsFrom(settings);
    	m_bwtIndex.loadSettingsFrom(settings);
    	m_readType.loadSettingsFrom(settings);
    	m_checkColorSpaced.loadSettingsFrom(settings);
    	m_checkIndexRefSeq.loadSettingsFrom(settings);
    	m_readGroup.loadSettingsFrom(settings);
    	m_readGroupBoolean.loadSettingsFrom(settings);
    	m_alnalgo.loadSettingsFrom(settings);
    	m_ALN_THREADS.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	/** added for HTE **/
    	super.validateSettings(settings);
    	
    	m_bwafile.validateSettings(settings);
    	m_refseqfile.validateSettings(settings);
    	m_bwtIndex.validateSettings(settings);
    	m_readType.validateSettings(settings);
    	m_checkColorSpaced.validateSettings(settings);
    	m_checkIndexRefSeq.validateSettings(settings);
    	m_readGroup.validateSettings(settings);
    	m_readGroupBoolean.validateSettings(settings);
    	m_alnalgo.validateSettings(settings);
    	m_ALN_THREADS.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
    	
    	
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {

    }
}


package de.helmholtz_muenchen.ibis.ngs.gatkrealignment;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.lofs.PathProcessor;

/**
 * This is the model implementation of GATKRealignment.
 * 
 *
 * @author
 */
public class GATKRealignmentNodeModel extends HTExecutorNodeModel {

	protected static final NodeLogger logger = NodeLogger.getLogger(GATKRealignmentNodeModel.class);

	// general options

	// path to gatk executable
	static final String CFGKEY_GATK = "gatk";
	private final SettingsModelString m_gatk = new SettingsModelString(CFGKEY_GATK, "");
	// path to reference genome
	static final String CFGKEY_REF_GENOME = "ref_genome";
	private final SettingsModelString m_ref_genome = new SettingsModelString(CFGKEY_REF_GENOME, "");
	// use 1000 genomes phase 1 indel set
	static final String CFGKEY_USE_PHASE1_1000G = "use_phase1_1000G";
	static final boolean DEF_USE_PHASE1_1000G = true;
	private final SettingsModelBoolean m_use_phase1_1000G = new SettingsModelBoolean(CFGKEY_USE_PHASE1_1000G,
			DEF_USE_PHASE1_1000G);
	// path to 1000 genomes phase 1 indel set
	static final String CFGKEY_PHASE1_1000G_FILE = "phase1_1000G_file";
	static final String DEF_PHASE1_1000G_FILE = "";
	private final SettingsModelString m_phase1_1000G_file = new SettingsModelString(CFGKEY_PHASE1_1000G_FILE,
			DEF_PHASE1_1000G_FILE);
	// use mills and 1000 genomes gold standard indel set
	static final String CFGKEY_USE_MILLS_1000G = "use_mills_1000G";
	static final boolean DEF_USE_MILLS_1000G = true;
	private final SettingsModelBoolean m_use_mills_1000G = new SettingsModelBoolean(CFGKEY_USE_MILLS_1000G,
			DEF_USE_MILLS_1000G);
	// path to mills and 1000 genomes gold standard indel set
	static final String CFGKEY_MILLS_1000G_FILE = "mills_1000G_file";
	static final String DEF_MILLS_1000G_FILE = "";
	private final SettingsModelString m_mills_1000G_file = new SettingsModelString(CFGKEY_MILLS_1000G_FILE,
			DEF_MILLS_1000G_FILE);
	// use interval
	static final String CFGKEY_USE_INTERVAL = "use_interval";
	static final boolean DEF_USE_INTERVAL = false;
	private final SettingsModelBoolean m_use_interval = new SettingsModelBoolean(CFGKEY_USE_INTERVAL, DEF_USE_INTERVAL);
	// path to interval file
	static final String CFGKEY_INTERVAL_FILE = "interval_file";
	static final String DEF_INTERVAL_FILE = "";
	private final SettingsModelString m_interval_file = new SettingsModelString(CFGKEY_INTERVAL_FILE,
			DEF_INTERVAL_FILE);
	// number of threads for target creator
	static final String CFGKEY_NUM_THREADS = "num_threads";
	static final int DEF_NUM_THREADS = 1;
	static final int MIN_NUM_THREADS = 1;
	static final int MAX_NUM_THREADS = Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_num_threads = new SettingsModelIntegerBounded(CFGKEY_NUM_THREADS,
			DEF_NUM_THREADS, MIN_NUM_THREADS, MAX_NUM_THREADS);

	static final String CFGKEY_JAVAMEMORY = "gatkmemory";
	static final int DEF_NUM_JAVAMEMORY = 8;
	static final int MIN_NUM_JAVAMEMORY = 1;
	static final int MAX_NUM_JAVAMEMORY = Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_gatk_java_memory = new SettingsModelIntegerBounded(CFGKEY_JAVAMEMORY,
			DEF_NUM_JAVAMEMORY, MIN_NUM_JAVAMEMORY, MAX_NUM_JAVAMEMORY);

	// target creator

	// maximal interval length for realignment
	static final String CFGKEY_MAX_INTERVAL = "max_interval";
	static final int DEF_MAX_INTERVAL = 500;
	static final int MIN_MAX_INTERVAL = 1;
	static final int MAX_MAX_INTERVAL = Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_interval = new SettingsModelIntegerBounded(CFGKEY_MAX_INTERVAL,
			DEF_MAX_INTERVAL, MIN_MAX_INTERVAL, MAX_MAX_INTERVAL);
	// minimum reads for entropy calculation
	static final String CFGKEY_MIN_READS = "min_reads";
	static final int DEF_MIN_READS = 4;
	static final int MIN_MIN_READS = 1;
	static final int MAX_MIN_READS = Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_min_reads = new SettingsModelIntegerBounded(CFGKEY_MIN_READS,
			DEF_MIN_READS, MIN_MIN_READS, MAX_MIN_READS);
	// window size for SNP/high entropy clusters
	static final String CFGKEY_WINDOW = "window";
	static final int DEF_WINDOW = 10;
	static final int MIN_WINDOW = 1;
	static final int MAX_WINDOW = Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_window = new SettingsModelIntegerBounded(CFGKEY_WINDOW, DEF_WINDOW,
			MIN_WINDOW, MAX_WINDOW);
	// mismatch fraction (ungapped alignment)
	static final String CFGKEY_MISMATCH = "mismatch";
	static final double DEF_MISMATCH = 0.0;
	static final double MIN_MISMATCH = DEF_MISMATCH;
	static final double MAX_MISMATCH = 1.0;
	private final SettingsModelDoubleBounded m_mismatch = new SettingsModelDoubleBounded(CFGKEY_MISMATCH, DEF_MISMATCH,
			MIN_MISMATCH, MAX_MISMATCH);
	static final String CFGKEY_TC_OPT_FLAGS = "opt_flags";
	public final SettingsModelOptionalString m_tc_opt_flags = new SettingsModelOptionalString(CFGKEY_TC_OPT_FLAGS, "", false);

	// indel realignment

	// consensus determination model
	static final String CFGKEY_CONSENSUS_MODEL = "consensus_model";
	static final String[] VALUES_CONSENSUS_MODEL = { "USE_READS", "KNOWNS_ONLY", "USE_SW" };
	static final String DEF_CONSENSUS_MODEL = VALUES_CONSENSUS_MODEL[0];
	private final SettingsModelString m_consensus_model = new SettingsModelString(CFGKEY_CONSENSUS_MODEL,
			DEF_CONSENSUS_MODEL);
	// significance threshold
	static final String CFGKEY_LOD_THRESHOLD = "lod_threshold";
	static final double DEF_LOD_THRESHOLD = 5.0;
	static final double MIN_LOD_THRESHOLD = 0.0;
	static final double MAX_LOD_THRESHOLD = Double.MAX_VALUE;
	private final SettingsModelDoubleBounded m_lod_threshold = new SettingsModelDoubleBounded(CFGKEY_LOD_THRESHOLD,
			DEF_LOD_THRESHOLD, MIN_LOD_THRESHOLD, MAX_LOD_THRESHOLD);
	// entropy threshold
	static final String CFGKEY_ENTROPY = "entropy";
	static final double DEF_ENTROPY = 0.15;
	static final double MIN_ENTROPY = 0.0;
	static final double MAX_ENTROPY = 1.0;
	private final SettingsModelDoubleBounded m_entropy = new SettingsModelDoubleBounded(CFGKEY_ENTROPY, DEF_ENTROPY,
			MIN_ENTROPY, MAX_ENTROPY);
	// max # consensuses to try
	static final String CFGKEY_MAX_CONSENSUSES = "max_consensuses";
	static final int DEF_MAX_CONSENSUSES = 30;
	static final int MIN_MAX_CONSENSUSES = 1;
	static final int MAX_MAX_CONSENSUSES = Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_consensuses = new SettingsModelIntegerBounded(
			CFGKEY_MAX_CONSENSUSES, DEF_MAX_CONSENSUSES, MIN_MAX_CONSENSUSES, MAX_MAX_CONSENSUSES);
	// max insert size for realignment
	static final String CFGKEY_MAX_ISIZE = "max_isize";
	static final int DEF_MAX_ISIZE = 3000;
	static final int MIN_MAX_ISIZE = 0;
	static final int MAX_MAX_ISIZE = Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_isize = new SettingsModelIntegerBounded(CFGKEY_MAX_ISIZE,
			DEF_MAX_ISIZE, MIN_MAX_ISIZE, MAX_MAX_ISIZE);
	// max positional movement of read during realignment
	static final String CFGKEY_MAX_POS_MOVE = "max_pos_move";
	static final int DEF_MAX_POS_MOVE = 200;
	static final int MIN_MAX_POS_MOVE = 1;
	static final int MAX_MAX_POS_MOVE = Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_pos_move = new SettingsModelIntegerBounded(CFGKEY_MAX_POS_MOVE,
			DEF_MAX_POS_MOVE, MIN_MAX_POS_MOVE, MAX_MAX_POS_MOVE);
	// max # reads for consensus determination
	static final String CFGKEY_MAX_READS_CONS = "max_reads_cons";
	static final int DEF_MAX_READS_CONS = 120;
	static final int MIN_MAX_READS_CONS = 1;
	static final int MAX_MAX_READS_CONS = Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_reads_cons = new SettingsModelIntegerBounded(CFGKEY_MAX_READS_CONS,
			DEF_MAX_READS_CONS, MIN_MAX_READS_CONS, MAX_MAX_READS_CONS);
	// max # reads for realignment
	static final String CFGKEY_MAX_READS_REALIGN = "max_reads_realign";
	static final int DEF_MAX_READS_REALIGN = 20000;
	static final int MIN_MAX_READS_REALIGN = 1;
	static final int MAX_MAX_READS_REALIGN = Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_reads_realign = new SettingsModelIntegerBounded(
			CFGKEY_MAX_READS_REALIGN, DEF_MAX_READS_REALIGN, MIN_MAX_READS_REALIGN, MAX_MAX_READS_REALIGN);
	// do not output old cigar string ?
	static final String CFGKEY_ALIGNMENT_TAG = "alignment_tag";
	static final boolean DEF_ALIGNMENT_TAG = false;
	private final SettingsModelBoolean m_alignment_tag = new SettingsModelBoolean(CFGKEY_ALIGNMENT_TAG,
			DEF_ALIGNMENT_TAG);
	static final String CFGKEY_IR_OPT_FLAGS = "opt_flags";
	public final SettingsModelOptionalString m_ir_opt_flags = new SettingsModelOptionalString(CFGKEY_IR_OPT_FLAGS, "", false);

	// position of input files in input table
	private int posBam;
	public static final String OUT_COL1 = "Path2RealignedBAM";

	// Network/Proxy options
	// public static final String CFGKEY_USEPROXY="useproxy";
	// public static final String CFGKEY_PROXYHOST="proxyhost";
	// public static final String CFGKEY_PROXYPORT="proxyport";
	// public static final String CFGKEY_USEPROXYAUTH="useproxyauth";
	// public static final String CFGKEY_PROXYUSER="proxyuser";
	// public static final String CFGKEY_PROXYPASSWORD="proxypassword";

	// private final SettingsModelBoolean m_useproxy = new
	// SettingsModelBoolean(GATKRealignmentNodeModel.CFGKEY_USEPROXY, false);
	// private final SettingsModelString m_proxyhost = new SettingsModelString(
	// GATKRealignmentNodeModel.CFGKEY_PROXYHOST,"");
	// private final SettingsModelString m_proxyport = new SettingsModelString(
	// GATKRealignmentNodeModel.CFGKEY_PROXYPORT,"");
	// private final SettingsModelBoolean m_useproxyauth = new
	// SettingsModelBoolean(GATKRealignmentNodeModel.CFGKEY_USEPROXYAUTH,
	// false);
	// private final SettingsModelString m_proxyuser = new SettingsModelString(
	// GATKRealignmentNodeModel.CFGKEY_PROXYUSER,"");
	// private final SettingsModelString m_proxypassword = new
	// SettingsModelString(
	// GATKRealignmentNodeModel.CFGKEY_PROXYPASSWORD,"");

	/**
     * Constructor for the node model.
     */
    protected GATKRealignmentNodeModel() {
        super(1, 1);
        
    	addSetting(m_gatk);
        addSetting(m_ref_genome);
    	addSetting(m_use_phase1_1000G);
    	addSetting(m_phase1_1000G_file);
    	addSetting(m_use_mills_1000G);
    	addSetting(m_mills_1000G_file);
    	addSetting(m_use_interval);
    	addSetting(m_interval_file);
    	addSetting(m_num_threads);
    	addSetting(m_gatk_java_memory);
       
       //target creator
    	addSetting(m_max_interval);
    	addSetting(m_min_reads);
    	addSetting(m_mismatch);
    	addSetting(m_window);
    	addSetting(m_tc_opt_flags);
       
       //indel realignment
    	addSetting(m_consensus_model);
    	addSetting(m_lod_threshold);
    	addSetting(m_entropy);
    	addSetting(m_max_consensuses);
    	addSetting(m_max_isize);
    	addSetting(m_max_pos_move);
    	addSetting(m_max_reads_cons);
    	addSetting(m_max_reads_realign);
    	addSetting(m_alignment_tag);
    	addSetting(m_ir_opt_flags);
    	
        // file chooser for interval file is disabled from the beginning
        m_interval_file.setEnabled(false);
        
        //Proxy options
//    	m_proxyhost.setEnabled(false);
//    	m_proxyport.setEnabled(false);
//    	m_useproxyauth.setEnabled(false);
//    	m_proxyuser.setEnabled(false);
//    	m_proxypassword.setEnabled(false);
    }

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec)
			throws Exception {

		String inputfile = inData[0].iterator().next().getCell(posBam).toString();
		if(CompatibilityChecker.inputFileNotOk(inputfile)) {
			throw new InvalidSettingsException("No BAM file in input table or BAM file does not exist!");
		}

		// process path to input file -> location and base name of output file
		String index_file = IO.replaceFileExtension(inputfile, "bai");		

		// check bam file index
		if (Files.notExists(Paths.get(index_file)) && Files.notExists(Paths.get(inputfile+".bai"))) {
			throw new InvalidSettingsException("Missing BAM file index for "+inputfile);
		}

//		String outint = PathProcessor.createOutputFile(base, "intervals", "realigned");
		String outint = IO.replaceFileExtension(inputfile, "realigned.intervals");
//		String outbam = PathProcessor.createOutputFile(base, "bam", "realigned");
		String outbam = IO.replaceFileExtension(inputfile, "realigned.bam");

//		 Enable proxy if needed
//		 String proxyOptions = "";
//		 if(m_useproxy.getBooleanValue()){
//		
//		 proxyOptions += " -Dhttp.proxyHost=" + m_proxyhost.getStringValue();
//		 proxyOptions += " -Dhttp.proxyPort=" + m_proxyport.getStringValue();
//		
//		 if(m_useproxyauth.getBooleanValue()){
//		
//		 proxyOptions += " -Dhttp.proxyUser=" + m_proxyuser.getStringValue();
//		 proxyOptions += " -Dhttp.proxyPassword=" +
//		 m_proxypassword.getStringValue();
//		 }
//		
//		 proxyOptions += " ";
//		 }

		// run target creator
		this.targetcreator(exec, inputfile, outint);

		// run realignment
		this.realign(exec, outint, outbam, inputfile);

		/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, BAMCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(outbam)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
        return new BufferedDataTable[]{outTable};
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {

		posBam = CompatibilityChecker.getFirstIndexCellType(inSpecs[0], "BAMCell");
    	if(posBam==-1) {
    		throw new InvalidSettingsException("This node is not compatible with the precedent node as there is no BAM file in the input table!");
    	}
    	
    	if(CompatibilityChecker.inputFileNotOk(m_gatk.getStringValue())) {
    		throw new InvalidSettingsException("Set path to the GenomeAnalysisTK.jar!");
    	}
    	
    	//check reference file
    	String reffile =  m_ref_genome.getStringValue();
    	if(CompatibilityChecker.inputFileNotOk(m_ref_genome.getStringValue())) {
    		throw new InvalidSettingsException("Set reference genome!");
    	}

		if (!Files.exists(Paths.get(reffile + ".fai"))) {
			throw new InvalidSettingsException("Reference sequence index: " + reffile + ".fai does not exist!");
		}

		String refbase = PathProcessor.getBase(reffile);
		if (!Files.exists(Paths.get(refbase + ".dict"))) {
			throw new InvalidSettingsException("Reference sequence dictionary: " + refbase + ".dict does not exist!");
		}
		
		//check data sets
		String phase1 = m_phase1_1000G_file.getStringValue();
		if(m_use_phase1_1000G.getBooleanValue() && CompatibilityChecker.inputFileNotOk(phase1)) {
			throw new InvalidSettingsException("Set 1000G Indel data set!");
		}
		
		if (m_use_phase1_1000G.getBooleanValue() && !Files.exists(Paths.get(phase1 + ".idx"))) {
			throw new InvalidSettingsException("1000G Indel index file: " + phase1 + ".idx does not exist!");
		}
		
		String mills = m_mills_1000G_file.getStringValue();
		if(m_use_mills_1000G.getBooleanValue() && CompatibilityChecker.inputFileNotOk(mills)) {
			throw new InvalidSettingsException("Set Mills data set!");
		}
		
		if (m_use_mills_1000G.getBooleanValue() && !Files.exists(Paths.get(mills + ".idx"))) {
			throw new InvalidSettingsException("Mills index file: " + mills + ".idx does not exist!");
		}
		
		if(m_use_interval.getBooleanValue() && CompatibilityChecker.inputFileNotOk(m_interval_file.getStringValue())) {
			throw new InvalidSettingsException("Interval file not specified or does not exist!");
		}

		DataColumnSpec[] colspec = {new DataColumnSpecCreator(OUT_COL1, BAMCell.TYPE).createSpec()};
		return new DataTableSpec[] { new DataTableSpec(colspec) };
	}
	
	private void targetcreator(ExecutionContext exec, String inputbam, String outputint) throws Exception{
		
		String lockFile = outputint + SuccessfulRunChecker.LOCK_ENDING;
		
		//create command string
		
		// each data thread requires 2 GB of memory
		int threads = m_num_threads.getIntValue();
		int memory = m_gatk_java_memory.getIntValue()*threads;

		String cmd="java -jar -Xmx"+memory+"G " + m_gatk.getStringValue();
		cmd+=" -T RealignerTargetCreator";
		cmd+=" -nt "+threads;
		cmd+=" -R "+m_ref_genome.getStringValue();
		cmd+=" -I "+inputbam;
		cmd+=" -o "+outputint;
		
		if(m_use_phase1_1000G.getBooleanValue()){
			cmd+=" -known "+m_phase1_1000G_file.getStringValue();
		}
		
		if(m_use_mills_1000G.getBooleanValue()){
			cmd+=" -known "+m_mills_1000G_file.getStringValue();
		}
		
		if(m_use_interval.getBooleanValue()){
			cmd+=" -L "+m_interval_file.getStringValue();
		}
		
		cmd+=" -maxInterval "+m_max_interval.getIntValue();
		cmd+=" -minReads "+m_min_reads.getIntValue();
		cmd+=" -mismatch "+m_mismatch.getDoubleValue();
		cmd+=" -window "+m_window.getIntValue();

		if(m_tc_opt_flags.isActive()) {
			cmd+=" "+m_tc_opt_flags.getStringValue();
		}
		
		GATKRealignmentNodeModel.logger.info("Running GATK TargetCreator...");
		GATKRealignmentNodeModel.logger.info("Log files can be found in "+outputint+".stdOut and "+outputint+".stdErr");
		
		super.executeCommand(new String[]{cmd}, exec, new File(lockFile),outputint+".stdOut", outputint+".stdErr");
	}
	
	private void realign (ExecutionContext exec, String outputint, String outputbam, String inputbam) throws Exception{
    	
		String lockFile = outputbam + SuccessfulRunChecker.LOCK_ENDING;

		String cmd="java -jar -Xmx"+m_gatk_java_memory.getIntValue()+"G " + m_gatk.getStringValue();
		cmd+=" -T IndelRealigner";
		cmd+=" -R "+m_ref_genome.getStringValue();
		cmd+=" -I "+inputbam;
		cmd+=" -o "+outputbam;
		cmd+=" -targetIntervals "+outputint;
		
		if(m_use_phase1_1000G.getBooleanValue()){
			cmd+=" -known "+m_phase1_1000G_file.getStringValue();
		}
		
		if(m_use_mills_1000G.getBooleanValue()){
			cmd+=" -known "+m_mills_1000G_file.getStringValue();
		}
		
		if(m_use_interval.getBooleanValue()){
			cmd+=" -L "+m_interval_file.getStringValue();
		}
		
		cmd+=" -model "+m_consensus_model.getStringValue();
		cmd+=" -LOD "+m_lod_threshold.getDoubleValue();
		cmd+=" -entropy "+m_entropy.getDoubleValue();
		cmd+=" -maxConsensuses "+m_max_consensuses.getIntValue();
		cmd+=" -maxIsize "+m_max_isize.getIntValue();
		cmd+=" -maxPosMove "+m_max_pos_move.getIntValue();
		cmd+=" -greedy "+m_max_reads_cons.getIntValue();
		cmd+=" -maxReads "+m_max_reads_realign.getIntValue();
		
		if(m_alignment_tag.getBooleanValue()){
			cmd+=" -noTags ";
		}
		
		if(m_ir_opt_flags.isActive()) {
			cmd+= " " + m_ir_opt_flags.getStringValue();
		}
		
		GATKRealignmentNodeModel.logger.info("Running GATK IndelRealigner...");
		GATKRealignmentNodeModel.logger.info("Log files can be found in "+outputbam+".stdOut and "+outputbam+".stdErr");
		
		super.executeCommand(new String[]{cmd}, exec, new File(lockFile),outputbam+".stdOut", outputbam+".stdErr");
	}
}
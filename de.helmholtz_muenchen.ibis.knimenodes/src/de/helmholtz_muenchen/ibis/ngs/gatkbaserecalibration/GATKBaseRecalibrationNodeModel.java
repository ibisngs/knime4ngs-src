/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.ngs.gatkbaserecalibration;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.lofs.PathProcessor;


/**
 * This is the model implementation of GATKBaseRecalibration.
 * @author Marie-Sophie Friedl
 * @author Maximilian Hastreiter
 * @author Tim Jeske
 */
public class GATKBaseRecalibrationNodeModel extends HTExecutorNodeModel {
    
    protected static final NodeLogger logger = NodeLogger.getLogger(GATKBaseRecalibrationNodeModel.class);
    
    //general options
    
    // path to gatk executable
    static final String CFGKEY_GATK="gatk";
    private final SettingsModelString m_gatk = new SettingsModelString(CFGKEY_GATK, "");
    // path to reference genome
    static final String CFGKEY_REF_GENOME = "ref_genome";
    private final SettingsModelString m_ref_genome = new SettingsModelString(CFGKEY_REF_GENOME,"");
    // use 1000 genomes phase 1 indel set
    static final String CFGKEY_USE_PHASE1_1000G="use_phase1_1000G";
    static final boolean DEF_USE_PHASE1_1000G=true;
	private final SettingsModelBoolean m_use_phase1_1000G=new SettingsModelBoolean(CFGKEY_USE_PHASE1_1000G, DEF_USE_PHASE1_1000G);
	// path to 1000 genomes phase 1 indel set
	static final String CFGKEY_PHASE1_1000G_FILE="phase1_1000G_file";
	static final String DEF_PHASE1_1000G_FILE="";
	private final SettingsModelString m_phase1_1000G_file=new SettingsModelString(CFGKEY_PHASE1_1000G_FILE, DEF_PHASE1_1000G_FILE);
	// use mills and 1000 genomes gold standard indel set
	static final String CFGKEY_USE_MILLS_1000G="use_mills_1000G";
	static final boolean DEF_USE_MILLS_1000G=true;
	private final SettingsModelBoolean m_use_mills_1000G=new SettingsModelBoolean(CFGKEY_USE_MILLS_1000G, DEF_USE_MILLS_1000G);
	// path to mills and 1000 genomes gold standard indel set
	static final String CFGKEY_MILLS_1000G_FILE="mills_1000G_file";
	static final String DEF_MILLS_1000G_FILE="";
	private final SettingsModelString m_mills_1000G_file=new SettingsModelString(CFGKEY_MILLS_1000G_FILE, DEF_MILLS_1000G_FILE);
	// use dbsnp snp set
	static final String CFGKEY_USE_DBSNP="use_dbsnp";
	static final boolean DEF_USE_DBSNP=true;
	private final SettingsModelBoolean m_use_dbsnp=new SettingsModelBoolean(CFGKEY_USE_DBSNP, DEF_USE_DBSNP);
	// path to dbsnp snp set
	static final String CFGKEY_DBSNP_FILE="dbsnp_file";
	static final String DEF_DBSNP_FILE="";
	private final SettingsModelString m_dbsnp_file=new SettingsModelString(CFGKEY_DBSNP_FILE, DEF_DBSNP_FILE);
    // use interval
    static final String CFGKEY_USE_INTERVAL="use_interval";
    static final boolean DEF_USE_INTERVAL=false;
	private final SettingsModelBoolean m_use_interval = new SettingsModelBoolean(CFGKEY_USE_INTERVAL, DEF_USE_INTERVAL);
    // path to interval file
    static final String CFGKEY_INTERVAL_FILE="interval_file";
    static final String DEF_INTERVAL_FILE="";
	private final SettingsModelString m_interval_file = new SettingsModelString(CFGKEY_INTERVAL_FILE, DEF_INTERVAL_FILE);
	// create before after plots
	static final String CFGKEY_CREATE_PLOTS="create_plots";
	static final boolean DEF_CREATE_PLOTS=false;
	private final SettingsModelBoolean m_create_plots = new SettingsModelBoolean(CFGKEY_CREATE_PLOTS, DEF_CREATE_PLOTS);
	// number of cpu threads
	static final String CFGKEY_CPU_THREADS="cpu_threads";
	static final int DEF_CPU_THREADS=1;
	static final int MIN_CPU_THREADS=1;
	static final int MAX_CPU_THREADS=Integer.MAX_VALUE;
	final SettingsModelIntegerBounded m_cpu_threads= new SettingsModelIntegerBounded(CFGKEY_CPU_THREADS, DEF_CPU_THREADS, MIN_CPU_THREADS, MAX_CPU_THREADS);
	
	public static final String CFGKEY_JAVAMEMORY = "gatkmemory";
    public static final int DEF_NUM_JAVAMEMORY=8;
    public static final int MIN_NUM_JAVAMEMORY=1;
    public static final int MAX_NUM_JAVAMEMORY=Integer.MAX_VALUE;
    private final SettingsModelIntegerBounded m_GATK_JAVA_MEMORY = new SettingsModelIntegerBounded(CFGKEY_JAVAMEMORY, DEF_NUM_JAVAMEMORY, MIN_NUM_JAVAMEMORY, MAX_NUM_JAVAMEMORY);
	
    //BaseRecalibrator +PrintReads options
	
	//calculate context covariate
//	static final String CFGKEY_CONTEXT_COV= "context_cov";
//	static final boolean DEF_CONTEXT_COV=true;
//	private final SettingsModelBoolean m_context_cov = new SettingsModelBoolean(CFGKEY_CONTEXT_COV, DEF_CONTEXT_COV);
	//calculate cycle covariate
//	static final String CFGKEY_CYCLE_COV="cycle_cov";
//	static final boolean DEF_CYCLE_COV=true;
//	private final SettingsModelBoolean m_cycle_cov = new SettingsModelBoolean(CFGKEY_CYCLE_COV, DEF_CYCLE_COV);
	// calculate repeat length covariate
//	static final String CFGKEY_REP_LEN_COV="rep_len_cov";
//	static final boolean DEF_REP_LEN_COV=false;
//	private final SettingsModelBoolean m_rep_len_cov = new SettingsModelBoolean(CFGKEY_REP_LEN_COV, DEF_REP_LEN_COV);
	// calculate repeat unit covariate
//	static final String CFGKEY_REP_UNIT_COV="rep_unit_cov";
//	static final boolean DEF_REP_UNIT_COV=false;
//	private final SettingsModelBoolean m_rep_unit_cov= new SettingsModelBoolean(CFGKEY_REP_UNIT_COV, DEF_REP_UNIT_COV);
	// minimum base quality of bases in read tail
	static final String CFGKEY_LOW_QUAL_TAIL="low_qual_tail";
	static final int DEF_LOW_QUAL_TAIL=2;
	static final int MIN_LOW_QUAL_TAIL=0;
	static final int MAX_LOW_QUAL_TAIL=Byte.MAX_VALUE;
	private final SettingsModelIntegerBounded m_low_qual_tail= new SettingsModelIntegerBounded(CFGKEY_LOW_QUAL_TAIL, DEF_LOW_QUAL_TAIL, MIN_LOW_QUAL_TAIL, MAX_LOW_QUAL_TAIL);
	// gap open penalty
	static final String CFGKEY_GAP_OPEN="gap_open";
	static final double DEF_GAP_OPEN=40.0;
	static final double MIN_GAP_OPEN=0;
	static final double MAX_GAP_OPEN=Double.MAX_VALUE;
	private final SettingsModelDoubleBounded m_gap_open = new SettingsModelDoubleBounded(CFGKEY_GAP_OPEN, DEF_GAP_OPEN, MIN_GAP_OPEN, MAX_GAP_OPEN);
	// deletion default quality
	static final String CFGKEY_DELETION_DEF_QUAL="deletion_def_qual";
	static final int DEF_DELETION_DEF_QUAL=45;
	static final int MIN_DELETION_DEF_QUAL=-1;
	static final int MAX_DELETION_DEF_QUAL=Byte.MAX_VALUE;
	private final SettingsModelIntegerBounded m_deletion_def_qual = new SettingsModelIntegerBounded(CFGKEY_DELETION_DEF_QUAL, DEF_DELETION_DEF_QUAL, MIN_DELETION_DEF_QUAL, MAX_DELETION_DEF_QUAL);
	// insertion default quality
	static final String CFGKEY_INSERTION_DEF_QUAL="insertion_def_qual";
	static final int DEF_INSERTION_DEF_QUAL=45;
	static final int MIN_INSERTION_DEF_QUAL=-1;
	static final int MAX_INSERTION_DEF_QUAL=Byte.MAX_VALUE;
	private final SettingsModelIntegerBounded m_insertion_def_qual = new SettingsModelIntegerBounded(CFGKEY_INSERTION_DEF_QUAL, DEF_INSERTION_DEF_QUAL, MIN_INSERTION_DEF_QUAL, MAX_INSERTION_DEF_QUAL);
	// indel context size
	static final String CFGKEY_INDEL_CONTEXT_SIZE="indel_context_size";
	static final int DEF_INDEL_CONTEXT_SIZE=3;
	static final int MIN_INDEL_CONTEXT_SIZE=1;
	static final int MAX_INDEL_CONTEXT_SIZE=13;	
	private final SettingsModelIntegerBounded m_indel_context_size= new SettingsModelIntegerBounded(CFGKEY_INDEL_CONTEXT_SIZE, DEF_INDEL_CONTEXT_SIZE, MIN_INDEL_CONTEXT_SIZE, MAX_INDEL_CONTEXT_SIZE);
	// default quality mismatches
	static final String CFGKEY_MISMATCH_DEF_QUAL="mismatch_def_qual";
	static final int DEF_MISMATCH_DEF_QUAL=-1;
	static final int MIN_MISMATCH_DEF_QUAL=DEF_MISMATCH_DEF_QUAL;
	static final int MAX_MISMACTH_DEF_QUAL=Byte.MAX_VALUE;
	private final SettingsModelIntegerBounded m_mismatch_def_qual = new SettingsModelIntegerBounded(CFGKEY_MISMATCH_DEF_QUAL, DEF_MISMATCH_DEF_QUAL, MIN_MISMATCH_DEF_QUAL, MAX_MISMACTH_DEF_QUAL);
	// context size mismatches
	static final String CFGKEY_MISMATCH_CONTEXT_SIZE="mismatch_context_size";
	static final int DEF_MISMATCH_CONTEXT_SIZE=2;
	static final int MIN_MISMATCH_CONTEXT_SIZE=1;
	static final int MAX_MISMATCH_CONTEXT_SIZE=13;
	private final SettingsModelIntegerBounded m_mismatch_context_size = new SettingsModelIntegerBounded(CFGKEY_MISMATCH_CONTEXT_SIZE, DEF_MISMATCH_CONTEXT_SIZE, MIN_MISMATCH_CONTEXT_SIZE, MAX_MISMATCH_CONTEXT_SIZE);
	//maximum #cycles
	static final String CFGKEY_MAX_CYCLES="max_cycles";
	static final int DEF_MAX_CYCLES=500;
	static final int MIN_MAX_CYCLES=0;
	static final int MAX_MAX_CYCLES=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_cycles = new SettingsModelIntegerBounded(CFGKEY_MAX_CYCLES, DEF_MAX_CYCLES, MIN_MAX_CYCLES, MAX_MAX_CYCLES);
	// do not output additional tags
	static final String CFGKEY_SIMPLIFY_OUT="simplify_out";
	static final boolean DEF_SIMPLIY_OUT=false;
	private final SettingsModelBoolean m_simplify_out=new SettingsModelBoolean(CFGKEY_SIMPLIFY_OUT, DEF_SIMPLIY_OUT);
	
	static final String CFGKEY_BR_OPT_FLAGS = "br_opt_flags";
	private final SettingsModelOptionalString m_br_opt_flags = new SettingsModelOptionalString(CFGKEY_BR_OPT_FLAGS,"",false);
	static final String CFGKEY_AC_OPT_FLAGS = "ac_opt_flags";
	private final SettingsModelOptionalString m_ac_opt_flags = new SettingsModelOptionalString(CFGKEY_AC_OPT_FLAGS,"",false);
	static final String CFGKEY_PR_OPT_FLAGS = "pr_opt_flags";
	private final SettingsModelOptionalString m_pr_opt_flags = new SettingsModelOptionalString(CFGKEY_PR_OPT_FLAGS,"",false);
	
	//Network/Proxy options
//	public static final String CFGKEY_USEPROXY="useproxy";
//	public static final String CFGKEY_PROXYHOST="proxyhost";
//	public static final String CFGKEY_PROXYPORT="proxyport";
//	public static final String CFGKEY_USEPROXYAUTH="useproxyauth";
//	public static final String CFGKEY_PROXYUSER="proxyuser";
//	public static final String CFGKEY_PROXYPASSWORD="proxypassword";
	
//	private final SettingsModelBoolean m_useproxy = new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USEPROXY, false);
//	private final SettingsModelString m_proxyhost = new SettingsModelString(
//			GATKBaseRecalibrationNodeModel.CFGKEY_PROXYHOST,"");
//	private final SettingsModelString m_proxyport = new SettingsModelString(
//			GATKBaseRecalibrationNodeModel.CFGKEY_PROXYPORT,"");
//	private final SettingsModelBoolean m_useproxyauth = new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USEPROXYAUTH, false);
//	private final SettingsModelString m_proxyuser = new SettingsModelString(
//			GATKBaseRecalibrationNodeModel.CFGKEY_PROXYUSER,"");
//	private final SettingsModelString m_proxypassword = new SettingsModelString(
//			GATKBaseRecalibrationNodeModel.CFGKEY_PROXYPASSWORD,"");
 
	public static String OUT_COL1 = "Path2RecalibratedBAM";
	private int posBam;
	private String gatk_jar, ref_genome, phase1, mills, dbsnp;
	
    /**
     * Constructor for the node model.
     */
    protected GATKBaseRecalibrationNodeModel() {
    
        super(1, 1, 1);
        
    	addSetting(m_gatk);
    	addSetting(m_ref_genome);
    	addSetting(m_use_phase1_1000G);
    	addSetting(m_phase1_1000G_file);
    	addSetting(m_use_mills_1000G);
    	addSetting(m_mills_1000G_file);
    	addSetting(m_use_dbsnp);
    	addSetting(m_dbsnp_file);    	
    	addSetting(m_use_interval);
    	addSetting(m_interval_file);
    	addSetting(m_create_plots);
    	addSetting(m_cpu_threads);
    	addSetting(m_GATK_JAVA_MEMORY);
    	
//    	addSetting(m_context_cov);
//    	addSetting(m_cycle_cov);
//    	addSetting(m_rep_len_cov);
//    	addSetting(m_rep_unit_cov);
    	addSetting(m_gap_open);
    	addSetting(m_low_qual_tail);
    	addSetting(m_max_cycles);
    	addSetting(m_indel_context_size);
    	addSetting(m_mismatch_context_size);
    	addSetting(m_deletion_def_qual);
    	addSetting(m_insertion_def_qual);
    	addSetting(m_mismatch_def_qual);
    	addSetting(m_simplify_out);
        addSetting(m_br_opt_flags);
        addSetting(m_ac_opt_flags);
        addSetting(m_pr_opt_flags);
    	
        m_interval_file.setEnabled(false);
        
        //Proxy options
//    	m_proxyhost.setEnabled(false);
//    	m_proxyport.setEnabled(false);
//    	m_useproxyauth.setEnabled(false);
//    	m_proxyuser.setEnabled(false);
//    	m_proxypassword.setEnabled(false);
        
		addPrefPageSetting(m_gatk, IBISKNIMENodesPlugin.GATK);
		addPrefPageSetting(m_ref_genome, IBISKNIMENodesPlugin.REF_GENOME);
		addPrefPageSetting(m_phase1_1000G_file, IBISKNIMENodesPlugin.RES_1000G_INDELS);
		addPrefPageSetting(m_mills_1000G_file, IBISKNIMENodesPlugin.RES_MILLS);
		addPrefPageSetting(m_dbsnp_file, IBISKNIMENodesPlugin.RES_DBSNP);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	//Check input table integrity
    	CompatibilityChecker.inDataCheck(inData);
    	
        // check bam input file
        String inputfile=inData[0].iterator().next().getCell(posBam).toString();
        
        if(CompatibilityChecker.inputFileNotOk(inputfile)) {
			throw new InvalidSettingsException("No BAM file in input table or BAM file does not exist!");
		}

		// process path to input file -> location and base name of output file
		String index_file = IO.replaceFileExtension(inputfile, "bai");

		// check bam file index
		if (Files.notExists(Paths.get(index_file)) && Files.notExists(Paths.get(inputfile+".bai"))) {
			throw new InvalidSettingsException("Missing BAM file index for "+inputfile);
		}
       
//        boolean [] covariates ={true, true, true, true, false, false};
//        if(!m_create_plots.getBooleanValue()){
//            covariates[0]=m_context_cov.getBooleanValue();
//            covariates[1]=m_cycle_cov.getBooleanValue();
//            covariates[2]=true;
//            covariates[3]=true;
//            covariates[4]=m_rep_len_cov.getBooleanValue();
//            covariates[5]=m_rep_unit_cov.getBooleanValue();
//        }
        
        int[] indelmis = new int[5];
        indelmis[0]= m_deletion_def_qual.getIntValue();
        indelmis[1]= m_insertion_def_qual.getIntValue();
        indelmis[2]= m_mismatch_def_qual.getIntValue();
		indelmis[3]= m_indel_context_size.getIntValue();
        indelmis[4]= m_mismatch_context_size.getIntValue();
        
        // file names for tool output
        String recaltable =IO.replaceFileExtension(inputfile, "recal.table");
        String recalbam = IO.replaceFileExtension(inputfile, "recal.bam");
        
		//Enable proxy if needed
//		String proxyOptions = "";
//		if(m_useproxy.getBooleanValue()){
//			
//			proxyOptions += " -Dhttp.proxyHost=" + m_proxyhost.getStringValue();
//			proxyOptions += " -Dhttp.proxyPort=" + m_proxyport.getStringValue();
//			
//			if(m_useproxyauth.getBooleanValue()){
//				
//    			proxyOptions += " -Dhttp.proxyUser=" + m_proxyuser.getStringValue();
//    			proxyOptions += " -Dhttp.proxyPassword=" + m_proxypassword.getStringValue();
//			}
//			
//			proxyOptions += " ";
//		}
		
    	runBaseRecalibrator(exec, inputfile, recaltable, indelmis);
    	
    	if(m_create_plots.getBooleanValue()){
    		
    		// additional output
    		String recalaftertable = IO.replaceFileExtension(inputfile, "post_recal.table");
    		String recalplots = IO.replaceFileExtension(inputfile, "recal_plots.pdf");
    		String recalintermediate = IO.replaceFileExtension(inputfile, "recal_plots_intermediateData.csv");
    		
    		runBaseRecalibrator(exec, inputfile, recaltable, recalaftertable, indelmis);
    		runAnalyzeCovariates(exec, recaltable, recalaftertable, recalplots, recalintermediate);
    	}
    	
    	runPrintReads(exec, inputfile, recaltable, recalbam);
    	
    	// create column specifications
    	DataColumnSpec [] colspec = {new DataColumnSpecCreator(OUT_COL1, BAMCell.TYPE).createSpec()};
	    DataTableSpec outspec=new DataTableSpec(colspec);
	    BufferedDataContainer c = exec.createDataContainer(outspec);
	    FileCell [] row = {FileCellFactory.create(recalbam)};
	    c.addRowToTable(new DefaultRow("row0", row));
	    c.close();
	    BufferedDataTable out=c.getTable();
 
        return new BufferedDataTable [] {out};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	super.updatePrefs();
    	gatk_jar = IO.processFilePath(m_gatk.getStringValue());
    	ref_genome = IO.processFilePath(m_ref_genome.getStringValue());
    	phase1 = IO.processFilePath(m_phase1_1000G_file.getStringValue());
    	mills = IO.processFilePath(m_mills_1000G_file.getStringValue());
    	dbsnp = IO.processFilePath(m_dbsnp_file.getStringValue());
    	
    	if(GATKBaseRecalibrationNodeDialog.getUseMainInputColBool()){
    		posBam = inSpecs[0].findColumnIndex(GATKBaseRecalibrationNodeDialog.getMainInputCol1());
    		if(!inSpecs[0].getColumnSpec(posBam).getType().toString().equals("BAMCell")){
    			posBam = -1;
    		}
    	} else {
    		posBam = CompatibilityChecker.getFirstIndexCellType(inSpecs[0], "BAMCell");
    	}
    	if(!(posBam>-1)) {
    		throw new InvalidSettingsException("This node is not compatible with the precedent node as there is no BAM file in the input table!");
    	}
    	
    	if(CompatibilityChecker.inputFileNotOk(gatk_jar)) {
    		throw new InvalidSettingsException("Set path to the GenomeAnalysisTK.jar!");
    	}
    	
    	//check reference file
    	if(CompatibilityChecker.inputFileNotOk(ref_genome)) {
    		throw new InvalidSettingsException("Set reference genome!");
    	}

		if (!Files.exists(Paths.get(ref_genome + ".fai"))) {
			throw new InvalidSettingsException("Reference sequence index: " + ref_genome + ".fai does not exist!");
		}

		String refbase = PathProcessor.getBase(ref_genome);
		if (!Files.exists(Paths.get(refbase + ".dict"))) {
			throw new InvalidSettingsException("Reference sequence dictionary: " + refbase + ".dict does not exist!");
		}
		
		//check data sets
		if(m_use_phase1_1000G.getBooleanValue() && CompatibilityChecker.inputFileNotOk(phase1)) {
			throw new InvalidSettingsException("Set 1000G Indel data set!");
		}
		
//		if (m_use_phase1_1000G.getBooleanValue() && !Files.exists(Paths.get(phase1 + ".idx"))) {
//			throw new InvalidSettingsException("1000G Indel index file: " + phase1 + ".idx does not exist!");
//		}
		
		if(m_use_mills_1000G.getBooleanValue() && CompatibilityChecker.inputFileNotOk(mills)) {
			throw new InvalidSettingsException("Set Mills data set!");
		}
		
//		if (m_use_mills_1000G.getBooleanValue() && !Files.exists(Paths.get(mills + ".idx"))) {
//			throw new InvalidSettingsException("Mills index file: " + mills + ".idx does not exist!");
//		}
		
		if(m_use_dbsnp.getBooleanValue() && CompatibilityChecker.inputFileNotOk(dbsnp)) {
			throw new InvalidSettingsException("Set dbSNP data set!");
		}
		
//		if (m_use_dbsnp.getBooleanValue() && !Files.exists(Paths.get(dbsnp + ".idx"))) {
//			throw new InvalidSettingsException("dbSNP index file: " + dbsnp + ".idx does not exist!");
//		}
		
		if(m_use_interval.getBooleanValue() && CompatibilityChecker.inputFileNotOk(IO.processFilePath(m_interval_file.getStringValue()))) {
			throw new InvalidSettingsException("Interval file not specified or does not exist!");
		}

		DataColumnSpec[] colspec = {new DataColumnSpecCreator(OUT_COL1, BAMCell.TYPE).createSpec()};
		return new DataTableSpec[] { new DataTableSpec(colspec) };
    }
    
    protected void runBaseRecalibrator(ExecutionContext exec, String inputbam, String outtable, int [] indelmis) throws Exception {
    	runBaseRecalibrator(exec, inputbam, null, outtable, indelmis);
    }
    
	protected void runBaseRecalibrator(ExecutionContext exec, String inputbam, String inputtable, String outtable, int [] indelmis) throws Exception {
		
		String lockFile = outtable + SuccessfulRunChecker.LOCK_ENDING;
		
		//create command string
		String cmd="java -jar -Xmx"+m_GATK_JAVA_MEMORY.getIntValue()+"G "+ gatk_jar;
		cmd+=" -T BaseRecalibrator";
		cmd+=" -R "+ref_genome;
		cmd+=" -I "+inputbam;
		if(inputtable!=null) {
			cmd+=" -BQSR "+inputtable;
		}
		cmd+=" -o "+outtable;
		
		if(m_use_phase1_1000G.getBooleanValue()){
			cmd+=" -knownSites "+phase1;
		}
		
		if(m_use_mills_1000G.getBooleanValue()){
			cmd+=" -knownSites "+mills;
		}
		
		if(m_use_dbsnp.getBooleanValue()){
			cmd+=" -knownSites "+dbsnp;
		}
		
		if(m_use_interval.getBooleanValue()){
			cmd+=" -L "+IO.processFilePath(m_interval_file.getStringValue());
		}
		
//		cmd+=" -noStandard";
//		if(cov[0]){
//			cmd+=" -cov ContextCovariate";
//		}
//		if(cov[1]){
//			cmd+=" -cov CycleCovariate";
//		}
//		if(cov[2]){
//			cmd+=" -cov QualityScoreCovariate";
//		}
//		if(cov[3]){
//			cmd+=" -cov ReadGroupCovariate";
//		}
//		if(cov[4]){
//			cmd+=" -cov RepeatLengthCovariate";
//		}
//		if(cov[5]){
//			cmd+=" -cov RepeatUnitCovariate";
//		}
		
		cmd+=" -lqt "+m_low_qual_tail.getIntValue();
		cmd+=" -bqsrBAQGOP "+m_gap_open.getDoubleValue();
		cmd+=" -maxCycle "+m_max_cycles.getIntValue();
		
		cmd+=" -ddq "+indelmis[0];
		cmd+=" -idq "+indelmis[1];
		cmd+=" -mdq "+indelmis[2];
		cmd+=" -ics "+indelmis[3];
		cmd+=" -mcs "+indelmis[4];
		
		cmd+=" -nct "+m_cpu_threads.getIntValue();
		
		if(m_br_opt_flags.isActive()) {
			cmd+=" "+m_br_opt_flags.getStringValue();
		}
		
		// run command
		GATKBaseRecalibrationNodeModel.logger.info("Running GATK BaseRecalibrator...");
		GATKBaseRecalibrationNodeModel.logger.info("Log files can be found in "+outtable+".stdOut and "+outtable+".stdErr");
		super.executeCommand(new String[]{cmd},outtable, exec, new File(lockFile),outtable+".stdOut", outtable+".stdErr");
	}
	
	protected void runPrintReads(ExecutionContext exec, String inputbam, String inputtable, String outbam) throws Exception {
		
		String lockFile = outbam + SuccessfulRunChecker.LOCK_ENDING;
		
		//create command string
		String cmd="java -jar -Xmx"+m_GATK_JAVA_MEMORY.getIntValue()+"G " + gatk_jar;
		cmd+=" -T PrintReads";
		cmd+=" -R "+ref_genome;
		cmd+=" -I "+inputbam;
		cmd+=" -BQSR "+inputtable;
		cmd+=" -o "+outbam;
		
		cmd+=" -nct "+m_cpu_threads.getIntValue();
		
		if(m_simplify_out.getBooleanValue()){
			cmd+=" -s ";
		}
		
		if(m_pr_opt_flags.isActive()) {
			cmd+=" " +m_pr_opt_flags.getStringValue();
		}
		
		GATKBaseRecalibrationNodeModel.logger.info("Running GATK PrintReads...");
		GATKBaseRecalibrationNodeModel.logger.info("Log files can be found in "+outbam+".stdOut and "+outbam+".stdErr");
		super.executeCommand(new String[]{cmd}, outbam, exec, new File(lockFile),outbam+".stdOut", outbam+".stdErr");
	}
	
	protected void runAnalyzeCovariates(ExecutionContext exec, String beforetable, String aftertable, String pplots, String recalintermediate) throws Exception {
		
		String lockFile = pplots + SuccessfulRunChecker.LOCK_ENDING;
		
		String cmd="java -jar -Xmx"+m_GATK_JAVA_MEMORY.getIntValue()+"G " + gatk_jar;
		cmd+=" -T AnalyzeCovariates";
		cmd+=" -R "+ref_genome;
		cmd+=" -before "+beforetable;
		cmd+=" -after "+aftertable;
		cmd+=" -plots "+pplots;
		cmd+=" -csv "+recalintermediate;
		
		if(m_use_interval.getBooleanValue()){
			cmd+=" -L "+IO.processFilePath(m_interval_file.getStringValue());
		}
		
		if(m_ac_opt_flags.isActive()) {
			cmd+=" " +m_ac_opt_flags.getStringValue();
		}
		
		// PATH environment is needed for calling Rscript
		Map<String, String> map =System.getenv();
		String [] env = new String []{"PATH="+map.get("PATH")};
		
				GATKBaseRecalibrationNodeModel.logger.info("Running GATK AnalyzeCovariates...");
		GATKBaseRecalibrationNodeModel.logger.info("Log files can be found in "+pplots+".stdOut and "+pplots+".stdErr");
		super.executeCommand(new String[]{cmd},  pplots, exec, env, new File(lockFile), pplots+".stdOut", pplots+".stdErr", null, null, null);
	}
}
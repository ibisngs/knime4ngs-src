package de.helmholtz_muenchen.ibis.ngs.gatkbaserecalibration;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;


import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.lofs.PathProcessor;


/**
 * This is the model implementation of GATKBaseRecalibration.
 * 
 *
 * @author 
 */
public class GATKBaseRecalibrationNodeModel extends HTExecutorNodeModel {
    
    protected static final NodeLogger logger = NodeLogger.getLogger(GATKBaseRecalibrationNodeModel.class);
    
    //general options
    
    // path to gatk executable
    static final String CFGKEY_GATK="gatk";
    private final SettingsModelString m_gatk = new SettingsModelString(CFGKEY_GATK, "");
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
	
	//BaseRecalibrator +PrintReads options
	
	//calculate context covariate
	static final String CFGKEY_CONTEXT_COV= "context_cov";
	static final boolean DEF_CONTEXT_COV=true;
	private final SettingsModelBoolean m_context_cov = new SettingsModelBoolean(CFGKEY_CONTEXT_COV, DEF_CONTEXT_COV);
	//calculate cycle covariate
	static final String CFGKEY_CYCLE_COV="cycle_cov";
	static final boolean DEF_CYCLE_COV=true;
	private final SettingsModelBoolean m_cycle_cov = new SettingsModelBoolean(CFGKEY_CYCLE_COV, DEF_CYCLE_COV);
	// calculate repeat length covariate
	static final String CFGKEY_REP_LEN_COV="rep_len_cov";
	static final boolean DEF_REP_LEN_COV=false;
	private final SettingsModelBoolean m_rep_len_cov = new SettingsModelBoolean(CFGKEY_REP_LEN_COV, DEF_REP_LEN_COV);
	// calculate repeat unit covariate
	static final String CFGKEY_REP_UNIT_COV="rep_unit_cov";
	static final boolean DEF_REP_UNIT_COV=false;
	private final SettingsModelBoolean m_rep_unit_cov= new SettingsModelBoolean(CFGKEY_REP_UNIT_COV, DEF_REP_UNIT_COV);
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
	
	
	static final String CFGKEY_OPT_FLAGS = "opt_flags";
	public final SettingsModelOptionalString m_opt_flags = new SettingsModelOptionalString(CFGKEY_OPT_FLAGS,"",false);
	
	private int posBam;
	private int posRef;
	
	//variables are only used when previous node is a gatk node
	boolean gatk=false;
	boolean p1=false;
	boolean mills=false;
	boolean dbsnp=false;
	private int posGatk;
	private int posP1;
	private int posMills;
	private int posDbsnp;

	//Network/Proxy options
	public static final String CFGKEY_USEPROXY="useproxy";
	public static final String CFGKEY_PROXYHOST="proxyhost";
	public static final String CFGKEY_PROXYPORT="proxyport";
	public static final String CFGKEY_USEPROXYAUTH="useproxyauth";
	public static final String CFGKEY_PROXYUSER="proxyuser";
	public static final String CFGKEY_PROXYPASSWORD="proxypassword";
	
	private final SettingsModelBoolean m_useproxy = new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USEPROXY, false);
	private final SettingsModelString m_proxyhost = new SettingsModelString(
			GATKBaseRecalibrationNodeModel.CFGKEY_PROXYHOST,"");
	private final SettingsModelString m_proxyport = new SettingsModelString(
			GATKBaseRecalibrationNodeModel.CFGKEY_PROXYPORT,"");
	private final SettingsModelBoolean m_useproxyauth = new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USEPROXYAUTH, false);
	private final SettingsModelString m_proxyuser = new SettingsModelString(
			GATKBaseRecalibrationNodeModel.CFGKEY_PROXYUSER,"");
	private final SettingsModelString m_proxypassword = new SettingsModelString(
			GATKBaseRecalibrationNodeModel.CFGKEY_PROXYPASSWORD,"");
	
	
	public static final String CFGKEY_JAVAMEMORY = "gatkmemory";
    public static final int DEF_NUM_JAVAMEMORY=8;
    public static final int MIN_NUM_JAVAMEMORY=1;
    public static final int MAX_NUM_JAVAMEMORY=Integer.MAX_VALUE;
    private final SettingsModelIntegerBounded m_GATK_JAVA_MEMORY = new SettingsModelIntegerBounded(CFGKEY_JAVAMEMORY, DEF_NUM_JAVAMEMORY, MIN_NUM_JAVAMEMORY, MAX_NUM_JAVAMEMORY);
 
	
	
    /**
     * Constructor for the node model.
     */
    protected GATKBaseRecalibrationNodeModel() {
    
        super(1, 1);
        
        m_interval_file.setEnabled(false);
        
        //Proxy options
    	m_proxyhost.setEnabled(false);
    	m_proxyport.setEnabled(false);
    	m_useproxyauth.setEnabled(false);
    	m_proxyuser.setEnabled(false);
    	m_proxypassword.setEnabled(false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

       	/*
       	 * input table:
    	 * check path to bam
    	 * bam format
    	 * index .bai available
    	 * check path to reference
    	 * reference.fa.fai
    	 * reference.dict
    	 * 
    	 * optional: gatk path
    	 * optional: phase1 indels
    	 * optional: mills and 1000g gold standard
    	 * optional: dbsnp
    	 */
    	
    	//retrieve information from table
        DataRow r=inData[0].iterator().next();
        
        // bam file
        
        // check bam input file
        String inputfile=r.getCell(posBam).toString();
        
        // check if path is null
        if(inputfile.equals("")){
        	throw new Exception("No bam file available, something went wrong with the previous node!");
        }
        
        // check path to bam file
        if(!Files.exists(Paths.get(inputfile))){
        	throw new Exception("Path to input bam file: "+inputfile+" does not exist");
        }
        
        //process path to input file -> location and base name of output file
        String base = PathProcessor.getBase(inputfile);
        String fileextension = PathProcessor.getExt(inputfile);

        // check bam format
        if(!fileextension.equals("bam")){
        	throw new Exception("Input file is not in bam format!");
        }
        
        // check bam file index
        if(!Files.exists(Paths.get(base+".bai"))){
        		throw new Exception("Missing bam file index: "+base+".bai");
        }
        
        // reference file
        
        String reffile=r.getCell(posRef).toString();
        
        // path to reffile should not be null
        if(reffile.equals("")){
        	throw new Exception("No reference file available, something went wrong with the previous node!");
        }
        
        // check path to reference path
        if(!Files.exists(Paths.get(reffile))){
        	throw new Exception("Reference sequence file: "+reffile+" does not exist");
        }
        
        // check path to reference index
        if(!Files.exists(Paths.get(reffile+".fai"))){
        	throw new Exception("Reference sequence index: "+reffile+".fai does not exist");
        }
        
        //process path to reference file
        String refbase=PathProcessor.getBase(reffile);
        
        // check path to reference sequence dictionary
        if(!Files.exists(Paths.get(refbase+".dict"))){
        	throw new Exception("Reference seuqnece dictionary: "+refbase+".dict does not exist");
        }
        
        //gatk executable
        
        String gatkfile="";
        // info from previous node
        if(gatk){
        	gatkfile=r.getCell(posGatk).toString();
            if(gatkfile.equals("")){
            	throw new Exception("No gatk executable available, something went wrong with the previous node!");
            }
        }
        // info from node dialog
        else{
        	gatkfile=m_gatk.getStringValue();
            if(gatkfile.equals("")){
            	throw new Exception("Missing GATK executable: You have to configure the node before executing it!");
            }  
        }
        
        if(!Files.exists(Paths.get(gatkfile))){
        	throw new Exception("GATK executable: "+gatkfile+" does not exist");
        }

        // phase 1 indels
                
        String phase1file="";
        if(m_use_phase1_1000G.getBooleanValue()){
        	// from previous node
        	if(p1){
        		phase1file=r.getCell(posP1).toString();
                if(phase1file.equals("")){
                	throw new Exception("No file with 1000 genomes phase1 indels available, something went wrong with the previous node!");
                }
        	}
        	// from node dialog
        	else{
        		phase1file=m_phase1_1000G_file.getStringValue();
                if(phase1file.equals("")){
                	throw new Exception("Missing 1000 genomes phase1 indel set: You have to configure the node before executing it!");
                } 
        	}
        	
        	//check path and index
            if(!Files.exists(Paths.get(phase1file))){
            	throw new Exception("1000 genomes phase1 indel file: "+phase1file+" does not exist");
            }
            if(!Files.exists(Paths.get(phase1file+".idx"))){
            	throw new Exception("1000 genomes phase1 indel index file: "+phase1file+".idx does not exist");
            }
        }
        

        
        // mills and 1000G gold standard
        
        String millsfile="";
        if(m_use_mills_1000G.getBooleanValue()){
        	// from previous node
        	if(mills){
        		millsfile=r.getCell(posMills).toString();
                if(millsfile.equals("")){
                	throw new Exception("No file with mills and 1000 genomes gold standard indels available, something went wrong with the previous node!");
                }
        	}
        	// from node dialog
        	else{
        		millsfile=m_mills_1000G_file.getStringValue();
                if(millsfile.equals("")){
                	throw new Exception("Missing mills and 1000 genomes gold standard indel set: You have to configure the node before executing it!");
                } 
        	}
        	
        	//check path and index
            if(!Files.exists(Paths.get(millsfile))){
            	throw new Exception("Mills and 1000 genomes gold standard indel file: "+millsfile+" does not exist");
            }
            if(!Files.exists(Paths.get(millsfile+".idx"))){
            	throw new Exception("Mills and 1000 genomes gold standard indel index file: "+millsfile+".idx does not exist");
            }
        }
        

        
        //dbsnp snps
    	
        String dbsnpfile="";
        if(m_use_dbsnp.getBooleanValue()){
        	if(dbsnp){
        		dbsnpfile=r.getCell(posDbsnp).toString();
                if(dbsnpfile.equals("")){
                	throw new Exception("No file with dbsnp snps available, something went wrong with the previous node!");
                }
        	}
        	else{
        		dbsnpfile=m_dbsnp_file.getStringValue();
                if(dbsnpfile.equals("")){
                	throw new Exception("Missing dbsnp snp set: You have to configure the node before executing it!");
                } 
        	}
        	
        	// check path and index
            if(!Files.exists(Paths.get(dbsnpfile))){
            	throw new Exception("Dbsnp file: "+dbsnpfile+" does not exist");
            }
            if(!Files.exists(Paths.get(dbsnpfile+".idx"))){
            	throw new Exception("Dbsnp index file: "+dbsnpfile+".idx does not exist");
            }
        }

        
        // interval file
        
        String intfile="";
        if(m_use_interval.getBooleanValue()){
        	intfile=m_interval_file.getStringValue();
            if(intfile.equals("")){
            	throw new Exception("Missing interval file: You have to configure the node properly!");
            }
            //check path
            if(!Files.exists(Paths.get(intfile))){
            	throw new Exception("Interval file: "+intfile+" does not exist");
            }
        }
       
        boolean [] covariates ={true, true, true, true, false, false};
        if(!m_create_plots.getBooleanValue()){
            covariates[0]=m_context_cov.getBooleanValue();
            covariates[1]=m_cycle_cov.getBooleanValue();
            covariates[2]=true;
            covariates[3]=true;
            covariates[4]=m_rep_len_cov.getBooleanValue();
            covariates[5]=m_rep_unit_cov.getBooleanValue();
        }
        
        int[] indelmis = new int[5];
        indelmis[0]= m_deletion_def_qual.getIntValue();
        indelmis[1]= m_insertion_def_qual.getIntValue();
        indelmis[2]= m_mismatch_def_qual.getIntValue();
		indelmis[3]= m_indel_context_size.getIntValue();
        indelmis[4]= m_mismatch_context_size.getIntValue();
        
        // file names for tool output
        String recaltable=PathProcessor.createOutputFile(base, "table", "recal");
        String recalbam= PathProcessor.createOutputFile(base, "bam", "recal");
        
		//Enable proxy if needed
		String proxyOptions = "";
		if(m_useproxy.getBooleanValue()){
			
			proxyOptions += " -Dhttp.proxyHost=" + m_proxyhost.getStringValue();
			proxyOptions += " -Dhttp.proxyPort=" + m_proxyport.getStringValue();
			
			if(m_useproxyauth.getBooleanValue()){
				
    			proxyOptions += " -Dhttp.proxyUser=" + m_proxyuser.getStringValue();
    			proxyOptions += " -Dhttp.proxyPassword=" + m_proxypassword.getStringValue();
			}
			
			proxyOptions += " ";
		}
        
		int GATK_MEMORY_USAGE = m_GATK_JAVA_MEMORY.getIntValue();
		
    	new RunGATKBaseRecalibration().BaseRecalibrator(exec, gatkfile, inputfile, reffile, recaltable, phase1file, millsfile, dbsnpfile, intfile, covariates, m_low_qual_tail.getIntValue(), m_gap_open.getDoubleValue(), m_max_cycles.getIntValue(), indelmis, m_cpu_threads.getIntValue(), proxyOptions, GATK_MEMORY_USAGE, m_opt_flags.getStringValue());
    	
    	if(m_create_plots.getBooleanValue()){
    		
    		// additional output
    		String recalaftertable=PathProcessor.createOutputFile(base, "table", "post_recal");
    		String recalplots=PathProcessor.createOutputFile(base, "pdf", "recal_plots");
    		String recalintermediate=PathProcessor.createOutputFile(base, "csv", "recal_plots_intermediateData");
    		
    		new RunGATKBaseRecalibration().BaseRecalibrator(exec, gatkfile, inputfile, reffile, recaltable, recalaftertable, phase1file, millsfile, dbsnpfile, intfile, covariates,m_low_qual_tail.getIntValue(), m_gap_open.getDoubleValue(), m_max_cycles.getIntValue(), indelmis, m_cpu_threads.getIntValue(), proxyOptions, GATK_MEMORY_USAGE, m_opt_flags.getStringValue());
    		new RunGATKBaseRecalibration().AnalyzeCovariates(exec, gatkfile, reffile, recaltable, recalaftertable, recalplots, intfile, proxyOptions, GATK_MEMORY_USAGE,recalintermediate);
    		
    	}
    	
    	new RunGATKBaseRecalibration().PrintReads(exec, gatkfile, inputfile, reffile, recaltable, recalbam, m_simplify_out.getBooleanValue(), m_cpu_threads.getIntValue(), proxyOptions, GATK_MEMORY_USAGE);
    	
    	/*
    	 * output table
    	 * path to bam
    	 * path to refseq
    	 * path to gatk
    	 * if set path to 100G phase1 indels
    	 * if set path to mills and 1000G gold standard indels
    	 * if set path to dbsnp
    	 */
    	
    	//determine number of output columns
    	int colcount=3;
    	if(p1 || m_use_phase1_1000G.getBooleanValue()){
    		colcount++;
    	}
    	if(mills || m_use_mills_1000G.getBooleanValue()){
    		colcount++;
    	}
    	if(dbsnp || m_use_dbsnp.getBooleanValue()){
    		colcount++;
    	}
    	
    	// create column specifications
    	DataColumnSpec [] colspec = new DataColumnSpec[colcount];
    	int count=0;
    	colspec[count++]=new DataColumnSpecCreator("Path2BAMFile", StringCell.TYPE).createSpec();
    	colspec[count++]=new DataColumnSpecCreator("Path2SEQFile", StringCell.TYPE).createSpec();
    	colspec[count++]=new DataColumnSpecCreator("Path2GATKFile", StringCell.TYPE).createSpec();
    	if(p1 || m_use_phase1_1000G.getBooleanValue()){
        	colspec[count++]=new DataColumnSpecCreator("Path2phase1", StringCell.TYPE).createSpec();
    	}
    	if(mills || m_use_mills_1000G.getBooleanValue()){
    		colspec[count++]=new DataColumnSpecCreator("Path2mills", StringCell.TYPE).createSpec();
    	}
    	if(dbsnp || m_use_dbsnp.getBooleanValue()){
    		colspec[count++]=new DataColumnSpecCreator("Path2dbsnp", StringCell.TYPE).createSpec();
    	}
    	
    	
    	//create table
	    DataTableSpec outspec=new DataTableSpec(colspec);
	    BufferedDataContainer c = exec.createDataContainer(outspec);
	    
	    // fill string cells
	    DataCell [] row = new DataCell [colcount];
	    count=0;
	    row[count++]=new StringCell(recalbam);
	    row[count++]=new StringCell(reffile);
	    row[count++]=new StringCell(gatkfile);
	    if(p1){
		    row[count++]=new StringCell(r.getCell(posP1).toString());
	    }
	    else if (m_use_phase1_1000G.getBooleanValue()){
	    	row[count++]=new StringCell(phase1file);
	    }
	    if(mills){
		    row[count++]=new StringCell(r.getCell(posMills).toString());
	    }
	    else if (m_use_mills_1000G.getBooleanValue()){
	    	row[count++]=new StringCell(millsfile);
	    }
    	if(dbsnp){
		    row[count++]=new StringCell(r.getCell(posDbsnp).toString());
    	}
    	else if (m_use_dbsnp.getBooleanValue()){
    		row[count++]=new StringCell(dbsnpfile);
    	}

	    //create row and add it to the container
	    c.addRowToTable(new DefaultRow("row0", row));
	    
	    //create final table
	    c.close();
	    BufferedDataTable out=c.getTable();
 
        return new BufferedDataTable [] {out};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
            	    	
    	//input port: BAMLoader, PicardTools, GATKRealignment
    	//reference file: "Sequence file", "Path2SEQFile", "Path2RefFile"
    	//input sam/bam file: "Path2BAMFile"
    	
    	String [] cols= inSpecs[0].getColumnNames();
    	
    	
    	//check if path to bam file is available
    	if(!inSpecs[0].containsName("Path2BAMFile")){
    		if(inSpecs[0].containsName("Path2SAMFile")){
    			throw new InvalidSettingsException("Your input file is in sam format! GATK requires an indexed bam file which is sorted by genomic coordinate. Try out the PicardTools node.");
    		}
    		else{
    			throw new InvalidSettingsException("Previous node is incompatible! Missing path to bam file!");
    		}
    	}
    	
    	//check if path to reference sequence file is available
    	if(!inSpecs[0].containsName("Path2SEQFile")){
    		throw new InvalidSettingsException("Previous node is incompatible! Missing path to reference sequence!");
    	}
    	
    	// if previous node is a gatk node -> pass informations/files on to this node
    	if(inSpecs[0].containsName("Path2GATKFile")){
    		logger.info("Previous node is a gatk node");
    		
    		gatk=true;
    		m_gatk.setEnabled(false);
    		
    		if(inSpecs[0].containsName("Path2phase1")){
    			p1=true;
    			m_phase1_1000G_file.setEnabled(false);
    		}
    		
    		if(inSpecs[0].containsName("Path2mills")){
    			mills=true;
    			m_mills_1000G_file.setEnabled(false);
    		}

    		if(inSpecs[0].containsName("Path2dbsnp")){
    			dbsnp=true;
    			m_dbsnp_file.setEnabled(false);
    		}
    	}
    	
    	for (int i=0; i<cols.length; i++){
    		if(cols[i].equals("Path2BAMFile")){
    			posBam=i;
    		}
    		if(cols[i].equals("Path2SEQFile")){
    			posRef=i;
    		}
    		if(cols[i].equals("Path2GATKFile")){
    			posGatk=i;
    		}
    		if(cols[i].equals("Path2phase1")){
    			posP1=i;
    		}
    		if(cols[i].equals("Path2mills")){
    			posMills=i;
    		}
    		if(cols[i].equals("Path2dbsnp")){
    			posDbsnp=i;
    		}
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

        // general options
    	m_gatk.saveSettingsTo(settings);
    	m_use_phase1_1000G.saveSettingsTo(settings);
    	m_phase1_1000G_file.saveSettingsTo(settings);
    	m_use_mills_1000G.saveSettingsTo(settings);
    	m_mills_1000G_file.saveSettingsTo(settings);
    	m_use_dbsnp.saveSettingsTo(settings);
    	m_dbsnp_file.saveSettingsTo(settings);
    	m_use_interval.saveSettingsTo(settings);
    	m_interval_file.saveSettingsTo(settings);
    	m_create_plots.saveSettingsTo(settings);
    	m_cpu_threads.saveSettingsTo(settings);
    	m_GATK_JAVA_MEMORY.saveSettingsTo(settings);
    	
    	m_context_cov.saveSettingsTo(settings);
    	m_cycle_cov.saveSettingsTo(settings);
    	m_rep_len_cov.saveSettingsTo(settings);
    	m_rep_unit_cov.saveSettingsTo(settings);
    	m_gap_open.saveSettingsTo(settings);
    	m_low_qual_tail.saveSettingsTo(settings);
    	m_max_cycles.saveSettingsTo(settings);
    	m_mismatch_context_size.saveSettingsTo(settings);
    	m_indel_context_size.saveSettingsTo(settings);
    	m_deletion_def_qual.saveSettingsTo(settings);
    	m_insertion_def_qual.saveSettingsTo(settings);
    	m_mismatch_def_qual.saveSettingsTo(settings);
    	m_simplify_out.saveSettingsTo(settings);
    	
    	//Proxy options
    	m_useproxy.saveSettingsTo(settings);
    	m_proxyhost.saveSettingsTo(settings);
    	m_proxyport.saveSettingsTo(settings);
    	m_useproxyauth.saveSettingsTo(settings);
    	m_proxyuser.saveSettingsTo(settings);
    	m_proxypassword.saveSettingsTo(settings);
    	
    	m_opt_flags.saveSettingsTo(settings);
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	/** added for HTE **/
    	super.loadValidatedSettingsFrom(settings);
            
    	m_gatk.loadSettingsFrom(settings);
    	m_use_phase1_1000G.loadSettingsFrom(settings);
    	m_phase1_1000G_file.loadSettingsFrom(settings);
    	m_use_mills_1000G.loadSettingsFrom(settings);
    	m_mills_1000G_file.loadSettingsFrom(settings);
    	m_use_dbsnp.loadSettingsFrom(settings);
    	m_dbsnp_file.loadSettingsFrom(settings);
    	m_use_interval.loadSettingsFrom(settings);
    	m_interval_file.loadSettingsFrom(settings);
    	m_create_plots.loadSettingsFrom(settings);
    	m_cpu_threads.loadSettingsFrom(settings);
    	m_GATK_JAVA_MEMORY.loadSettingsFrom(settings);
    	
    	m_context_cov.loadSettingsFrom(settings);
    	m_cycle_cov.loadSettingsFrom(settings);
    	m_rep_len_cov.loadSettingsFrom(settings);
    	m_rep_unit_cov.loadSettingsFrom(settings);
    	m_gap_open.loadSettingsFrom(settings);
    	m_low_qual_tail.loadSettingsFrom(settings);
    	m_max_cycles.loadSettingsFrom(settings);
    	m_mismatch_context_size.loadSettingsFrom(settings);
    	m_indel_context_size.loadSettingsFrom(settings);
    	m_deletion_def_qual.loadSettingsFrom(settings);
    	m_insertion_def_qual.loadSettingsFrom(settings);
    	m_mismatch_def_qual.loadSettingsFrom(settings);
    	m_simplify_out.loadSettingsFrom(settings);
    	
    	//Proxy options
    	m_useproxy.loadSettingsFrom(settings);
    	m_proxyhost.loadSettingsFrom(settings);
    	m_proxyport.loadSettingsFrom(settings);
    	m_useproxyauth.loadSettingsFrom(settings);
    	m_proxyuser.loadSettingsFrom(settings);
    	m_proxypassword.loadSettingsFrom(settings);

    	m_opt_flags.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    
    	/** added for HTE **/
    	super.validateSettings(settings);
    	
    	m_gatk.validateSettings(settings);
    	m_use_phase1_1000G.validateSettings(settings);
    	m_phase1_1000G_file.validateSettings(settings);
    	m_use_mills_1000G.validateSettings(settings);
    	m_mills_1000G_file.validateSettings(settings);
    	m_use_dbsnp.validateSettings(settings);
    	m_dbsnp_file.validateSettings(settings);    	
    	m_use_interval.validateSettings(settings);
    	m_interval_file.validateSettings(settings);
    	m_create_plots.validateSettings(settings);
    	m_cpu_threads.validateSettings(settings);
    	m_GATK_JAVA_MEMORY.validateSettings(settings);
    	
    	m_context_cov.validateSettings(settings);
    	m_cycle_cov.validateSettings(settings);
    	m_rep_len_cov.validateSettings(settings);
    	m_rep_unit_cov.validateSettings(settings);
    	m_gap_open.validateSettings(settings);
    	m_low_qual_tail.validateSettings(settings);
    	m_max_cycles.validateSettings(settings);
    	m_indel_context_size.validateSettings(settings);
    	m_mismatch_context_size.validateSettings(settings);
    	m_deletion_def_qual.validateSettings(settings);
    	m_insertion_def_qual.validateSettings(settings);
    	m_mismatch_def_qual.validateSettings(settings);
    	m_simplify_out.validateSettings(settings);
    	
    	//Proxy options
    	m_useproxy.validateSettings(settings);
    	m_proxyhost.validateSettings(settings);
    	m_proxyport.validateSettings(settings);
    	m_useproxyauth.validateSettings(settings);
    	m_proxyuser.validateSettings(settings);
    	m_proxypassword.validateSettings(settings);
    	
    	m_opt_flags.validateSettings(settings);
    	
    	// check if at least one set of polymorphisms is used
    	if(!settings.getBoolean(CFGKEY_USE_PHASE1_1000G) && !settings.getBoolean(CFGKEY_USE_MILLS_1000G) && !settings.getBoolean(CFGKEY_USE_DBSNP)){
    		throw new InvalidSettingsException("You have to use at least one file containing polymorphisms");
    	}

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


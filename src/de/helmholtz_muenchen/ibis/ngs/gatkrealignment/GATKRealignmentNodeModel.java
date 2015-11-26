package de.helmholtz_muenchen.ibis.ngs.gatkrealignment;

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
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
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
    // use interval
    static final String CFGKEY_USE_INTERVAL="use_interval";
    static final boolean DEF_USE_INTERVAL=false;
	private final SettingsModelBoolean m_use_interval = new SettingsModelBoolean(CFGKEY_USE_INTERVAL, DEF_USE_INTERVAL);
    // path to interval file
    static final String CFGKEY_INTERVAL_FILE="interval_file";
    static final String DEF_INTERVAL_FILE="";
	private final SettingsModelString m_interval_file = new SettingsModelString(CFGKEY_INTERVAL_FILE, DEF_INTERVAL_FILE);
    // number of threads for target creator
    static final String CFGKEY_NUM_THREADS="num_threads";
    static final int DEF_NUM_THREADS=1;
    static final int MIN_NUM_THREADS=1;
    static final int MAX_NUM_THREADS=Integer.MAX_VALUE;
    private final SettingsModelIntegerBounded m_num_threads= new SettingsModelIntegerBounded(CFGKEY_NUM_THREADS, DEF_NUM_THREADS, MIN_NUM_THREADS, MAX_NUM_THREADS);
    
    static final String CFGKEY_JAVAMEMORY = "gatkmemory";
    static final int DEF_NUM_JAVAMEMORY=8;
    static final int MIN_NUM_JAVAMEMORY=1;
    static final int MAX_NUM_JAVAMEMORY=Integer.MAX_VALUE;
    private final SettingsModelIntegerBounded m_gatk_java_memory = new SettingsModelIntegerBounded(CFGKEY_JAVAMEMORY, DEF_NUM_JAVAMEMORY, MIN_NUM_JAVAMEMORY, MAX_NUM_JAVAMEMORY);
    
    
    // target creator
    
    // maximal interval length for realignment
    static final String CFGKEY_MAX_INTERVAL="max_interval";
    static final int DEF_MAX_INTERVAL=500;
    static final int MIN_MAX_INTERVAL=1;
    static final int MAX_MAX_INTERVAL=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_interval = new SettingsModelIntegerBounded(CFGKEY_MAX_INTERVAL, DEF_MAX_INTERVAL, MIN_MAX_INTERVAL, MAX_MAX_INTERVAL);
	// minimum reads for entropy calculation
	static final String CFGKEY_MIN_READS="min_reads";
	static final int DEF_MIN_READS=4;
	static final int MIN_MIN_READS=1;
	static final int MAX_MIN_READS=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_min_reads = new SettingsModelIntegerBounded(CFGKEY_MIN_READS, DEF_MIN_READS, MIN_MIN_READS, MAX_MIN_READS);
	// window size for SNP/high entropy clusters
	static final String CFGKEY_WINDOW="window";
	static final int DEF_WINDOW=10;
	static final int MIN_WINDOW=1;
	static final int MAX_WINDOW=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_window = new SettingsModelIntegerBounded(CFGKEY_WINDOW, DEF_WINDOW, MIN_WINDOW, MAX_WINDOW);
	// mismatch fraction (ungapped alignment)
	static final String CFGKEY_MISMATCH="mismatch";
	static final double DEF_MISMATCH=0.0;
	static final double MIN_MISMATCH=DEF_MISMATCH;
	static final double MAX_MISMATCH=1.0;
	private final SettingsModelDoubleBounded m_mismatch = new SettingsModelDoubleBounded(CFGKEY_MISMATCH, DEF_MISMATCH, MIN_MISMATCH, MAX_MISMATCH);
    
    // indel realignment
	
	//consensus determination model
	static final String CFGKEY_CONSENSUS_MODEL="consensus_model";
	static final String [] VALUES_CONSENSUS_MODEL={"USE_READS", "KNOWNS_ONLY", "USE_SW"};
	static final String DEF_CONSENSUS_MODEL=VALUES_CONSENSUS_MODEL[0];
	private final SettingsModelString m_consensus_model = new SettingsModelString(CFGKEY_CONSENSUS_MODEL, DEF_CONSENSUS_MODEL);
	//significance threshold
	static final String CFGKEY_LOD_THRESHOLD="lod_threshold";
	static final double DEF_LOD_THRESHOLD=5.0;
	static final double MIN_LOD_THRESHOLD=0.0;
	static final double MAX_LOD_THRESHOLD=Double.MAX_VALUE;
	private final SettingsModelDoubleBounded m_lod_threshold = new SettingsModelDoubleBounded(CFGKEY_LOD_THRESHOLD, DEF_LOD_THRESHOLD, MIN_LOD_THRESHOLD, MAX_LOD_THRESHOLD);
	// entropy threshold
	static final String CFGKEY_ENTROPY="entropy";
	static final double DEF_ENTROPY=0.15;
	static final double MIN_ENTROPY=0.0;
	static final double MAX_ENTROPY=1.0;
	private final SettingsModelDoubleBounded m_entropy = new SettingsModelDoubleBounded(CFGKEY_ENTROPY, DEF_ENTROPY, MIN_ENTROPY, MAX_ENTROPY);
	// max # consensuses to try
	static final String CFGKEY_MAX_CONSENSUSES="max_consensuses";
	static final int DEF_MAX_CONSENSUSES=30;
	static final int MIN_MAX_CONSENSUSES=1;
	static final int MAX_MAX_CONSENSUSES=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_consensuses = new SettingsModelIntegerBounded(CFGKEY_MAX_CONSENSUSES, DEF_MAX_CONSENSUSES, MIN_MAX_CONSENSUSES, MAX_MAX_CONSENSUSES);
	// max insert size for realignment
	static final String CFGKEY_MAX_ISIZE="max_isize";
	static final int DEF_MAX_ISIZE=3000;
	static final int MIN_MAX_ISIZE=0;
	static final int MAX_MAX_ISIZE=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_isize = new SettingsModelIntegerBounded(CFGKEY_MAX_ISIZE, DEF_MAX_ISIZE, MIN_MAX_ISIZE, MAX_MAX_ISIZE);
	// max positional movement of read during realignment
	static final String CFGKEY_MAX_POS_MOVE="max_pos_move";
	static final int DEF_MAX_POS_MOVE=200;
	static final int MIN_MAX_POS_MOVE=1;
	static final int MAX_MAX_POS_MOVE=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_pos_move = new SettingsModelIntegerBounded(CFGKEY_MAX_POS_MOVE, DEF_MAX_POS_MOVE, MIN_MAX_POS_MOVE, MAX_MAX_POS_MOVE);
	// max # reads for consensus determination
	static final String CFGKEY_MAX_READS_CONS="max_reads_cons";
	static final int DEF_MAX_READS_CONS=120;
	static final int MIN_MAX_READS_CONS=1;
	static final int MAX_MAX_READS_CONS=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_reads_cons = new SettingsModelIntegerBounded(CFGKEY_MAX_READS_CONS, DEF_MAX_READS_CONS, MIN_MAX_READS_CONS, MAX_MAX_READS_CONS);
	// max # reads for realignment
	static final String CFGKEY_MAX_READS_REALIGN="max_reads_realign";
	static final int DEF_MAX_READS_REALIGN=20000;
	static final int MIN_MAX_READS_REALIGN=1;
	static final int MAX_MAX_READS_REALIGN=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_reads_realign = new SettingsModelIntegerBounded(CFGKEY_MAX_READS_REALIGN, DEF_MAX_READS_REALIGN, MIN_MAX_READS_REALIGN, MAX_MAX_READS_REALIGN);
	// do not output old cigar string ?
	static final String CFGKEY_ALIGNMENT_TAG="alignment_tag";
	static final boolean DEF_ALIGNMENT_TAG=false;
	private final SettingsModelBoolean m_alignment_tag = new SettingsModelBoolean(CFGKEY_ALIGNMENT_TAG, DEF_ALIGNMENT_TAG);
    
    // position of input files in input table
    private int posBam;
    private int posRef;
    
	//variables are only used when previous node is a gatk node
	public static boolean gatk=false;
	public static boolean p1=false;
	public static boolean mills=false;
	public static boolean dbsnp=false;
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
	
	private final SettingsModelBoolean m_useproxy = new SettingsModelBoolean(GATKRealignmentNodeModel.CFGKEY_USEPROXY, false);
	private final SettingsModelString m_proxyhost = new SettingsModelString(
			GATKRealignmentNodeModel.CFGKEY_PROXYHOST,"");
	private final SettingsModelString m_proxyport = new SettingsModelString(
			GATKRealignmentNodeModel.CFGKEY_PROXYPORT,"");
	private final SettingsModelBoolean m_useproxyauth = new SettingsModelBoolean(GATKRealignmentNodeModel.CFGKEY_USEPROXYAUTH, false);
	private final SettingsModelString m_proxyuser = new SettingsModelString(
			GATKRealignmentNodeModel.CFGKEY_PROXYUSER,"");
	private final SettingsModelString m_proxypassword = new SettingsModelString(
			GATKRealignmentNodeModel.CFGKEY_PROXYPASSWORD,"");
		
	static final String CFGKEY_OPT_FLAGS = "opt_flags";
	public final SettingsModelOptionalString m_opt_flags = new SettingsModelOptionalString(CFGKEY_OPT_FLAGS,"",false);
	
    /**
     * Constructor for the node model.
     */
    protected GATKRealignmentNodeModel() {
    	
    	// 1 in port, 1 out port
        super(1, 1);
        

        // file chooser for interval file is disabled from the beginning
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
    protected BufferedDataTable[] execute (final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	
    	//retrieve information
        DataRow r=inData[0].iterator().next();
        
        
        //bam file
        
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
        
        // check bam file index
        if(!Files.exists(Paths.get(base+".bai"))){
        		throw new Exception("Missing bam file index: "+base+".bai");
        }
        
        // check bam format
        if(!fileextension.equals("bam")){
        	throw new Exception("Input file is not in bam format!");
        }
        
        
        // reference file
        String reffile = "";
        if(gatk) {
        	reffile=r.getCell(posRef).toString();
        } else {
        	reffile=m_ref_genome.getStringValue();
        }

        // path to reference should not be null
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
        String refbase = PathProcessor.getBase(reffile);
        
        // check path to reference sequence dictionary
        if(!Files.exists(Paths.get(refbase+".dict"))){
        	throw new Exception("Reference seuqnece dictionary: "+refbase+".dict does not exist");
        }
        
        
        //gatk executable
        
        String gatkfile;
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
        
        
        //phase1 indels
        
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
        
        
        // mills indels
    	
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
    	
    	
    	// create paths to output files
    	String outint =PathProcessor.createOutputFile(base, "intervals", "realigned");
        String outbam =PathProcessor.createOutputFile(base, "bam", "realigned");
        
        
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
        
		//GATK Memory Usage
		int GATK_MEMORY_USAGE = m_gatk_java_memory.getIntValue();
        
    	logger.info("Start GATK Realignment");
    	
    	// run target creator
    	 new RunGATKRealignment().targetcreator(exec, outint, inputfile, reffile, gatkfile, phase1file, millsfile, intfile, m_num_threads.getIntValue(), m_max_interval.getIntValue(), m_min_reads.getIntValue(), m_mismatch.getDoubleValue(), m_window.getIntValue(), proxyOptions, GATK_MEMORY_USAGE);
    	
    	//run realignment
    	 new RunGATKRealignment().realign(exec, outint, outbam, inputfile, reffile, gatkfile, phase1file, millsfile, intfile, m_consensus_model.getStringValue(), m_lod_threshold.getDoubleValue(), m_entropy.getDoubleValue(), m_max_consensuses.getIntValue(), m_max_isize.getIntValue(), m_max_pos_move.getIntValue(), m_max_reads_cons.getIntValue(), m_max_reads_realign.getIntValue(), m_alignment_tag.getBooleanValue(), proxyOptions, GATK_MEMORY_USAGE, m_opt_flags.getStringValue());
    	
    	// write output table
    	/*
    	 * output table
    	 * path to bam
    	 * path to refseq
    	 * path to gatk
    	 * if set path to 100G phase1 indels
    	 * if set path to mills and 1000G gold standard indels
    	 * if from previous node available path to dbsnp
    	 */
    	
    	//determine number of output columns
    	int colcount=3;
    	if(p1 || m_use_phase1_1000G.getBooleanValue()){
    		colcount++;
    	}
    	if(mills || m_use_mills_1000G.getBooleanValue()){
    		colcount++;
    	}
    	if(dbsnp){
    		colcount++;
    	}
    	
    	// create column specifications
    	DataColumnSpec [] colspec = new DataColumnSpec[colcount];
    	int count=0;
    	colspec[count++]=new DataColumnSpecCreator("Path2BAMFile", BAMCell.TYPE).createSpec();
    	colspec[count++]=new DataColumnSpecCreator("Path2SEQFile", FileCell.TYPE).createSpec();
    	colspec[count++]=new DataColumnSpecCreator("Path2GATKFile", FileCell.TYPE).createSpec();
    	if(p1 || m_use_phase1_1000G.getBooleanValue()){
        	colspec[count++]=new DataColumnSpecCreator("Path2phase1", FileCell.TYPE).createSpec();
    	}
    	if(mills || m_use_mills_1000G.getBooleanValue()){
    		colspec[count++]=new DataColumnSpecCreator("Path2mills", FileCell.TYPE).createSpec();
    	}
    	if(dbsnp){
    		colspec[count++]=new DataColumnSpecCreator("Path2dbsnp", FileCell.TYPE).createSpec();
    	}
    	
    	//create table
	    DataTableSpec outspec=new DataTableSpec(colspec);
	    BufferedDataContainer c = exec.createDataContainer(outspec);
	    
	    // fill string cells
	    DataCell [] row = new FileCell [colcount];
	    count=0;
	    row[count++]=(FileCell)FileCellFactory.create(outbam);
	    row[count++]=(FileCell)FileCellFactory.create(reffile);
	    row[count++]=(FileCell)FileCellFactory.create(gatkfile);
	    if(p1){
		    row[count++]=(FileCell)FileCellFactory.create(r.getCell(posP1).toString());
	    }
	    else if(m_use_phase1_1000G.getBooleanValue()){
	    	row[count++]=(FileCell)FileCellFactory.create(phase1file);
	    }
	    if(mills){
		    row[count++]=(FileCell)FileCellFactory.create(r.getCell(posMills).toString());
	    }
	    else if(m_use_mills_1000G.getBooleanValue()){
	    	row[count++]=(FileCell)FileCellFactory.create(millsfile);
	    }
    	if(dbsnp){
		    row[count++]=(FileCell)FileCellFactory.create(r.getCell(posDbsnp).toString());
    	}
	    
	    //create row and add it to the container
	    c.addRowToTable(new DefaultRow("row0", row));
	    
	    //create final table
	    c.close();
	    BufferedDataTable out=c.getTable();
	    
	    return new BufferedDataTable[]{out};

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
        
    	//input port: BAMLoader, PicardTools
    	//reference file: "Sequence file", "Path2SEQFile", "Path2RefFile"
    	//input sam/bam file: "Path2BAMFile"
    	
    	String [] cols= inSpecs[0].getColumnNames();
    	
    	//check if path to bam file is available
    	if(!CompatibilityChecker.checkInputCellType(inSpecs[0],"BAMCell")){
    		if(CompatibilityChecker.checkInputCellType(inSpecs[0],"SAMCell")){
    			throw new InvalidSettingsException("Your input file is in sam format! GATK requires an indexed bam file which is sorted by genomic coordinate. Try out the PicardTools node.");
    		}
    		else{
    			throw new InvalidSettingsException("Previous node is incompatible! Missing path to bam file!");
    		}
    	}
//    	
//    	//check if path to reference sequence file is available
//    	if(!inSpecs[0].containsName("Path2SEQFile")){
//    		throw new InvalidSettingsException("Previous node is incompatible! Missing path to reference sequence!");
//    	}
    	
    	// if previous node is a gatk node -> pass informations/files on to this node
    	if(inSpecs[0].containsName("Path2GATKFile")){
    		logger.info("Previous node is a gatk node");
    		
    		gatk=true;
    		m_gatk.setEnabled(false);
    		m_gatk.setStringValue("TableInput");
    		
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
    		}
    	}
    	
    	//get positions of input columns
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
        //general options
    	m_gatk.saveSettingsTo(settings);
        m_ref_genome.saveSettingsTo(settings);
    	m_use_phase1_1000G.saveSettingsTo(settings);
    	m_phase1_1000G_file.saveSettingsTo(settings);
    	m_use_mills_1000G.saveSettingsTo(settings);
    	m_mills_1000G_file.saveSettingsTo(settings);
    	m_use_interval.saveSettingsTo(settings);
    	m_interval_file.saveSettingsTo(settings);
    	m_num_threads.saveSettingsTo(settings);
    	m_gatk_java_memory.saveSettingsTo(settings);
       
       //target creator
    	m_max_interval.saveSettingsTo(settings);
    	m_min_reads.saveSettingsTo(settings);
    	m_mismatch.saveSettingsTo(settings);
    	m_window.saveSettingsTo(settings);
       
       //indel realignment
    	m_consensus_model.saveSettingsTo(settings);
    	m_lod_threshold.saveSettingsTo(settings);
    	m_entropy.saveSettingsTo(settings);
    	m_max_consensuses.saveSettingsTo(settings);
    	m_max_isize.saveSettingsTo(settings);
    	m_max_pos_move.saveSettingsTo(settings);
    	m_max_reads_cons.saveSettingsTo(settings);
    	m_max_reads_realign.saveSettingsTo(settings);
    	m_alignment_tag.saveSettingsTo(settings);
    	
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
    	// general options
    	m_gatk.loadSettingsFrom(settings);
        m_ref_genome.loadSettingsFrom(settings);
    	m_use_phase1_1000G.loadSettingsFrom(settings);
    	m_phase1_1000G_file.loadSettingsFrom(settings);
    	m_use_mills_1000G.loadSettingsFrom(settings);
    	m_mills_1000G_file.loadSettingsFrom(settings);
    	m_use_interval.loadSettingsFrom(settings);
    	m_interval_file.loadSettingsFrom(settings);
    	m_num_threads.loadSettingsFrom(settings);
    	m_gatk_java_memory.loadSettingsFrom(settings);
    	
    	// target creator
    	m_max_interval.loadSettingsFrom(settings);
    	m_min_reads.loadSettingsFrom(settings);
    	m_mismatch.loadSettingsFrom(settings);
    	m_window.loadSettingsFrom(settings);
    	
    	// indel realignment
    	m_consensus_model.loadSettingsFrom(settings);
    	m_lod_threshold.loadSettingsFrom(settings);
    	m_entropy.loadSettingsFrom(settings);
    	m_max_consensuses.loadSettingsFrom(settings);
    	m_max_isize.loadSettingsFrom(settings);
    	m_max_pos_move.loadSettingsFrom(settings);
    	m_max_reads_cons.loadSettingsFrom(settings);
    	m_max_reads_realign.loadSettingsFrom(settings);
    	m_alignment_tag.loadSettingsFrom(settings);
    	
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
    	//general options
    	
    	m_gatk.validateSettings(settings);
    	m_ref_genome.validateSettings(settings);
    	m_use_phase1_1000G.validateSettings(settings);
    	m_phase1_1000G_file.validateSettings(settings);    	
    	m_use_mills_1000G.validateSettings(settings);
    	m_mills_1000G_file.validateSettings(settings);    	
    	m_use_interval.validateSettings(settings);
    	m_interval_file.validateSettings(settings);
    	m_num_threads.validateSettings(settings);
    	m_gatk_java_memory.validateSettings(settings);
    	
    	//target creator
    	m_max_interval.validateSettings(settings);
    	m_min_reads.validateSettings(settings);
    	m_mismatch.validateSettings(settings);
    	m_window.validateSettings(settings);
    	
    	// indel realignment
    	m_consensus_model.validateSettings(settings);
    	m_lod_threshold.validateSettings(settings);
    	m_entropy.validateSettings(settings);
    	m_max_consensuses.validateSettings(settings);
    	m_max_isize.validateSettings(settings);
    	m_max_pos_move.validateSettings(settings);
    	m_max_reads_cons.validateSettings(settings);
    	m_max_reads_realign.validateSettings(settings);
    	m_alignment_tag.validateSettings(settings);
    	
    	//Proxy options
    	m_useproxy.validateSettings(settings);
    	m_proxyhost.validateSettings(settings);
    	m_proxyport.validateSettings(settings);
    	m_useproxyauth.validateSettings(settings);
    	m_proxyuser.validateSettings(settings);
    	m_proxypassword.validateSettings(settings);
    	
    	m_opt_flags.validateSettings(settings);
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
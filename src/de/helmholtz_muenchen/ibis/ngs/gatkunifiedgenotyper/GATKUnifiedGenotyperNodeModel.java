package de.helmholtz_muenchen.ibis.ngs.gatkunifiedgenotyper;

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
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.lofs.PathProcessor;


/**
 * This is the model implementation of GATKUnifiedGenotyper.
 * 
 *
 * @author 
 */
public class GATKUnifiedGenotyperNodeModel extends NodeModel {
    

    protected static final NodeLogger logger = NodeLogger.getLogger(GATKUnifiedGenotyperNodeModel.class);
    
    
    // path to gatk executable
    static final String CFGKEY_GATK="gatk";
    private final SettingsModelString m_gatk = new SettingsModelString(CFGKEY_GATK, "");
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
 	// call snps, indels or both
 	static final String CFGKEY_VARIANT_TYPE="variant_type";
 	static final String [] AVAIL_VARIANT_TYPE={"SNP", "INDEL", "SNP+INDEL"};
 	static final String DEF_VARIANT_TYPE=AVAIL_VARIANT_TYPE[2];
	private final SettingsModelString m_variant_type= new SettingsModelString(CFGKEY_VARIANT_TYPE, DEF_VARIANT_TYPE);
    // number of threads for target creator
    static final String CFGKEY_NUM_THREADS="num_threads";
    static final int DEF_NUM_THREADS=1;
    static final int MIN_NUM_THREADS=1;
    static final int MAX_NUM_THREADS=Integer.MAX_VALUE;
    private final SettingsModelIntegerBounded m_num_threads= new SettingsModelIntegerBounded(CFGKEY_NUM_THREADS, DEF_NUM_THREADS, MIN_NUM_THREADS, MAX_NUM_THREADS);
    
    // additional options
    
    // confidence threshold for calling
    static final String CFGKEY_CALL_MIN_CONFIDENCE="call_min_confidence";
    static final double DEF_CALL_MIN_CONFIDENCE=30.0;
    static final double MIN_CALL_MIN_CONFIDENCE=0.0;
    static final double MAX_CALL_MIN_CONFIDENCE=Double.MAX_VALUE;
	private final SettingsModelDoubleBounded m_call_min_confidence = new SettingsModelDoubleBounded(CFGKEY_CALL_MIN_CONFIDENCE, DEF_CALL_MIN_CONFIDENCE, MIN_CALL_MIN_CONFIDENCE, MAX_CALL_MIN_CONFIDENCE);
	// confidence threshold for emitting
	static final String CFGKEY_EMIT_MIN_CONFIDENCE="emit_min_confidence";
	static final double DEF_EMIT_MIN_CONFIDENCE=30.0;
    static final double MIN_EMIT_MIN_CONFIDENCE=0.0;
    static final double MAX_EMIT_MIN_CONFIDENCE=Double.MAX_VALUE;
	private final SettingsModelDoubleBounded m_emit_min_confidence = new SettingsModelDoubleBounded(CFGKEY_EMIT_MIN_CONFIDENCE, DEF_EMIT_MIN_CONFIDENCE, MIN_EMIT_MIN_CONFIDENCE, MAX_EMIT_MIN_CONFIDENCE);
    // pcr error rate
	static final String CFGKEY_PCR_ERR = "pcr_error";
	static final double DEF_PCR_ERR=1e-4;
	static final double MIN_PCR_ERR = 0;
	static final double MAX_PCR_ERR=1;
	private final SettingsModelDoubleBounded m_pcr_error = new SettingsModelDoubleBounded(CFGKEY_PCR_ERR, DEF_PCR_ERR, MIN_PCR_ERR, MAX_PCR_ERR);
	// fraction of contamination -> reads to filter
	static final String CFGKEY_CONTAMINATION = "contamination";
	static final double DEF_CONTAMINATION=0;
	static final double MIN_CONTAMINATION=0;
	static final double MAX_CONTAMINATION=1;
	private final SettingsModelDoubleBounded m_contamination = new SettingsModelDoubleBounded(CFGKEY_CONTAMINATION, DEF_CONTAMINATION, MIN_CONTAMINATION, MAX_CONTAMINATION);
	// heterozygosity value -> how much 2 individuals differ
	static final String CFGKEY_HET="heterozygosity";
	static final double DEF_HET=0.001;
	static final double MIN_HET=0;
	static final double MAX_HET=1;
	private final SettingsModelDoubleBounded m_heterozygosity = new SettingsModelDoubleBounded(CFGKEY_HET, DEF_HET, MIN_HET, MAX_HET);
	// fraction of deletions -> threshold for calling
	static final String CFGKEY_MAX_DELETION_FRAC="max_deletion_fraction";
	static final double DEF_MAX_DELETION_FRAC=0.05;
	static final double MIN_MAX_DELETION_FRAC=0;
	static final double MAX_MAX_DELETION_FRAC=1.1;
	private final SettingsModelDoubleBounded m_max_deletion_fraction = new SettingsModelDoubleBounded(CFGKEY_MAX_DELETION_FRAC, DEF_MAX_DELETION_FRAC, MIN_MAX_DELETION_FRAC, MAX_MAX_DELETION_FRAC);
	// min base quality threshold
	static final String CFGKEY_MIN_BASE_QUAL="min_base_qual";
	static final int DEF_MIN_BASE_QUAL=17;
	static final int MIN_MIN_BASE_QUAL=0;
	static final int MAX_MIN_BASE_QUAL=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_min_base_qual = new SettingsModelIntegerBounded(CFGKEY_MIN_BASE_QUAL, DEF_MIN_BASE_QUAL, MIN_MIN_BASE_QUAL, MAX_MIN_BASE_QUAL);
	// indel heterozygosity
	static final String CFGKEY_INDEL_HET="indel_heterozygosity";
	static final double DEF_INDEL_HET=1.25e-4;
	static final double MIN_INDEL_HET=0;
	static final double MAX_INDEL_HET=1;
	private final SettingsModelDoubleBounded m_indel_heterozygosity = new SettingsModelDoubleBounded(CFGKEY_INDEL_HET, DEF_INDEL_HET, MIN_INDEL_HET, MAX_INDEL_HET);
	// minimum indel count for calling indel
	static final String CFGKEY_MIN_INDEL_CNT="min_indel_count";
	static final int DEF_MIN_INDEL_CNT=5;
	static final int MIN_MIN_INDEL_CNT=0;
	static final int MAX_MIN_INDEL_CNT=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_min_indel_count = new SettingsModelIntegerBounded(CFGKEY_MIN_INDEL_CNT, DEF_MIN_INDEL_CNT, MIN_MIN_INDEL_CNT, MAX_MIN_INDEL_CNT);
	// minimum indel fraction for calling indel -> only sites fullfilling this criteria will be considered in count threshold
	static final String CFGKEY_MIN_INDEL_FRAC="min_indel_frac";
	static final double DEF_MIN_INDEL_FRAC=0.25;
	static final double MIN_MIN_INDEL_FRAC=0;
	static final double MAX_MIN_INDEL_FRAC=1;
	private final SettingsModelDoubleBounded m_min_indel_frac = new SettingsModelDoubleBounded(CFGKEY_MIN_INDEL_FRAC, DEF_MIN_INDEL_FRAC, MIN_MIN_INDEL_FRAC, MAX_MIN_INDEL_FRAC);
	// indel gap open penalty
	static final String CFGKEY_GAP_OPEN_PEN = "gap_open_pen";
	static final int DEF_GAP_OPEN_PEN= 45;
	static final int MIN_GAP_OPEN_PEN=0;
	static final int MAX_GAP_OPEN_PEN=Byte.MAX_VALUE;
	private final SettingsModelIntegerBounded m_gap_open_pen = new SettingsModelIntegerBounded(CFGKEY_GAP_OPEN_PEN, DEF_GAP_OPEN_PEN, MIN_GAP_OPEN_PEN, MAX_GAP_OPEN_PEN);
	// indel continuation penalty
	static final String CFGKEY_GAP_CONT_PEN="gap_cont_pen";
	static final int DEF_GAP_CONT_PEN=10;
	static final int MIN_GAP_CONT_PEN=0;
	static final int MAX_GAP_CONT_PEN=Byte.MAX_VALUE;
	private final SettingsModelIntegerBounded m_gap_cont_pen = new SettingsModelIntegerBounded(CFGKEY_GAP_CONT_PEN, DEF_GAP_CONT_PEN, MIN_GAP_CONT_PEN, MAX_GAP_CONT_PEN);
	// per-base alignment qualities
	static final String CFGKEY_BAQ="baq";
	static final String[] AVAIL_BAQ={"OFF", "CALCULATE_AS_NECESSARY", "RECALCULATE"};
	static final String DEF_BAQ=AVAIL_BAQ[1];
	private final SettingsModelString m_baq = new SettingsModelString(CFGKEY_BAQ, DEF_BAQ);
	// malformed read filter MBQ = mismatching base and quals
	static final String CFGKEY_MBQ="mbq";
	static final boolean DEF_MBQ=true;
	private final SettingsModelBoolean m_mbq = new SettingsModelBoolean(CFGKEY_MBQ, DEF_MBQ);
	
	private int posBam;
	private int posRef;
	
	//variables are only used when previous node is a gatk node
	boolean gatk=false;
	boolean p1=false;
	boolean mills=false;
	boolean dbsnp=false;
	private int posGatk;
	private int posMills;
	private int posP1;
	private int posDbsnp;
	
	//Network/Proxy options
	public static final String CFGKEY_USEPROXY="useproxy";
	public static final String CFGKEY_PROXYHOST="proxyhost";
	public static final String CFGKEY_PROXYPORT="proxyport";
	public static final String CFGKEY_USEPROXYAUTH="useproxyauth";
	public static final String CFGKEY_PROXYUSER="proxyuser";
	public static final String CFGKEY_PROXYPASSWORD="proxypassword";
	
	private final SettingsModelBoolean m_useproxy = new SettingsModelBoolean(GATKUnifiedGenotyperNodeModel.CFGKEY_USEPROXY, false);
	private final SettingsModelString m_proxyhost = new SettingsModelString(
			GATKUnifiedGenotyperNodeModel.CFGKEY_PROXYHOST,"");
	private final SettingsModelString m_proxyport = new SettingsModelString(
			GATKUnifiedGenotyperNodeModel.CFGKEY_PROXYPORT,"");
	private final SettingsModelBoolean m_useproxyauth = new SettingsModelBoolean(GATKUnifiedGenotyperNodeModel.CFGKEY_USEPROXYAUTH, false);
	private final SettingsModelString m_proxyuser = new SettingsModelString(
			GATKUnifiedGenotyperNodeModel.CFGKEY_PROXYUSER,"");
	private final SettingsModelString m_proxypassword = new SettingsModelString(
			GATKUnifiedGenotyperNodeModel.CFGKEY_PROXYPASSWORD,"");
	
    /**
     * Constructor for the node model.
     */
    protected GATKUnifiedGenotyperNodeModel() {
    
        // #in ports, #out ports
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
        String refbase = PathProcessor.getBase(reffile);
        
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
        //checks file path
        if(!Files.exists(Paths.get(gatkfile))){
        	throw new Exception("GATK executable: "+gatkfile+" does not exist");
        }

        //dbsnp snps
    	
        String dbsnpfile="";
        //dbsnp is used
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
        
        double [] param = new double [12];
        param[0]=m_call_min_confidence.getDoubleValue();
        param[1]=m_emit_min_confidence.getDoubleValue();
        param[2]=m_pcr_error.getDoubleValue();
        param[3]=m_contamination.getDoubleValue();
        param[4]=m_heterozygosity.getDoubleValue();
        param[5]=m_max_deletion_fraction.getDoubleValue();
        param[6]=m_min_base_qual.getIntValue();
        param[7]=m_indel_heterozygosity.getDoubleValue();
        param[8]=m_min_indel_count.getIntValue();
        param[9]=m_min_indel_frac.getDoubleValue();
        param[10]=m_gap_open_pen.getIntValue();
        param[11]=m_gap_cont_pen.getIntValue();
        
        
        // create file names of vcf output  
        String snpout="";
        if(m_variant_type.getStringValue().equals(AVAIL_VARIANT_TYPE[0]) || m_variant_type.getStringValue().equals(AVAIL_VARIANT_TYPE[2])){
        	snpout=PathProcessor.createOutputFile(base, "vcf", "snps");
        }
        
        String indelout="";
        if(m_variant_type.getStringValue().equals(AVAIL_VARIANT_TYPE[1]) || m_variant_type.getStringValue().equals(AVAIL_VARIANT_TYPE[2])){
        	indelout=PathProcessor.createOutputFile(base, "vcf", "indels");
        }
        
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
        
        
        // call snps
        if(!snpout.equals("")){
        	RunGATKUnifiedGenotyper.CallVariants(exec, gatkfile, inputfile, reffile, snpout, intfile, dbsnpfile, "SNP", m_num_threads.getIntValue(), param, m_baq.getStringValue(), m_mbq.getBooleanValue(), proxyOptions);
        }
        
        // call indels
        if(!indelout.equals("")){
        	RunGATKUnifiedGenotyper.CallVariants(exec, gatkfile, inputfile, reffile, indelout, intfile, dbsnpfile, "INDEL", m_num_threads.getIntValue(), param, m_baq.getStringValue(), m_mbq.getBooleanValue(), proxyOptions);
        }
        
    	//determine number of output columns
    	int colcount=4;
    	if(p1){
    		colcount++;
    	}
    	if(mills){
    		colcount++;
    	}
    	if(dbsnp || m_use_dbsnp.getBooleanValue()){
    		colcount++;
    	}
    	if(m_variant_type.getStringValue().equals(AVAIL_VARIANT_TYPE[2])){
    		colcount++;
    	}
    	
    	// create column specifications
    	DataColumnSpec [] colspec = new DataColumnSpec[colcount];
    	int count=0;
    	if(!snpout.equals("")){
        	colspec[count++]=new DataColumnSpecCreator("Path2VCFsnpFile", StringCell.TYPE).createSpec();
    	}
    	if(!indelout.equals("")){
        	colspec[count++]=new DataColumnSpecCreator("Path2VCFindelFile", StringCell.TYPE).createSpec();
    	}
    	colspec[count++]=new DataColumnSpecCreator("Path2BAMFile", StringCell.TYPE).createSpec();
    	colspec[count++]=new DataColumnSpecCreator("Path2SEQFile", StringCell.TYPE).createSpec();
    	colspec[count++]=new DataColumnSpecCreator("Path2GATKFile", StringCell.TYPE).createSpec();
    	if(p1){
        	colspec[count++]=new DataColumnSpecCreator("Path2phase1", StringCell.TYPE).createSpec();
    	}
    	if(mills){
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
    	if(!snpout.equals("")){
    	    row[count++]=new StringCell(snpout);
    	}
    	if(!indelout.equals("")){
    	    row[count++]=new StringCell(indelout);
    	}
	    row[count++]=new StringCell(inputfile);
	    row[count++]=new StringCell(reffile);
	    row[count++]=new StringCell(gatkfile);
	    if(p1){
		    row[count++]=new StringCell(r.getCell(posP1).toString());
	    }
	    if(mills){
		    row[count++]=new StringCell(r.getCell(posMills).toString());
	    }
    	if(dbsnp){
		    row[count++]=new StringCell(r.getCell(posDbsnp).toString());
    	}
    	else if(m_use_dbsnp.getBooleanValue()){
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
        // Code executed on reset.
        // Models build during execute are cleared here.
        // Also data handled in load/saveInternals will be erased here.
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
    		}
    		
    		if(inSpecs[0].containsName("Path2mills")){
    			mills=true;
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

        m_gatk.saveSettingsTo(settings);
        m_use_dbsnp.saveSettingsTo(settings);
        m_dbsnp_file.saveSettingsTo(settings);
        m_use_interval.saveSettingsTo(settings);
        m_interval_file.saveSettingsTo(settings);
        m_variant_type.saveSettingsTo(settings);
        m_num_threads.saveSettingsTo(settings);
        
        m_call_min_confidence.saveSettingsTo(settings);
        m_emit_min_confidence.saveSettingsTo(settings);
        m_pcr_error.saveSettingsTo(settings);
        m_contamination.saveSettingsTo(settings);
        m_heterozygosity.saveSettingsTo(settings);
        m_max_deletion_fraction.saveSettingsTo(settings);
        m_min_base_qual.saveSettingsTo(settings);
        m_indel_heterozygosity.saveSettingsTo(settings);
        m_min_indel_count.saveSettingsTo(settings);
        m_min_indel_frac.saveSettingsTo(settings);
        m_gap_open_pen.saveSettingsTo(settings);
        m_gap_cont_pen.saveSettingsTo(settings);
        m_baq.saveSettingsTo(settings);
        m_mbq.saveSettingsTo(settings);
        
    	//Proxy options
    	m_useproxy.saveSettingsTo(settings);
    	m_proxyhost.saveSettingsTo(settings);
    	m_proxyport.saveSettingsTo(settings);
    	m_useproxyauth.saveSettingsTo(settings);
    	m_proxyuser.saveSettingsTo(settings);
    	m_proxypassword.saveSettingsTo(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
    	m_gatk.loadSettingsFrom(settings);
    	m_use_dbsnp.loadSettingsFrom(settings);
    	m_dbsnp_file.loadSettingsFrom(settings);
    	m_use_interval.loadSettingsFrom(settings);
    	m_interval_file.loadSettingsFrom(settings);
    	m_variant_type.loadSettingsFrom(settings);
    	m_num_threads.loadSettingsFrom(settings);
    	
    	m_call_min_confidence.loadSettingsFrom(settings);
    	m_emit_min_confidence.loadSettingsFrom(settings);
    	m_pcr_error.loadSettingsFrom(settings);
    	m_contamination.loadSettingsFrom(settings);
    	m_heterozygosity.loadSettingsFrom(settings);
    	m_max_deletion_fraction.loadSettingsFrom(settings);
    	m_min_base_qual.loadSettingsFrom(settings);
    	m_indel_heterozygosity.loadSettingsFrom(settings);
    	m_min_indel_count.loadSettingsFrom(settings);
    	m_min_indel_frac.loadSettingsFrom(settings);
    	m_gap_open_pen.loadSettingsFrom(settings);
    	m_gap_cont_pen.loadSettingsFrom(settings);
    	m_baq.loadSettingsFrom(settings);
    	m_mbq.loadSettingsFrom(settings);
    	
    	//Proxy options
    	m_useproxy.loadSettingsFrom(settings);
    	m_proxyhost.loadSettingsFrom(settings);
    	m_proxyport.loadSettingsFrom(settings);
    	m_useproxyauth.loadSettingsFrom(settings);
    	m_proxyuser.loadSettingsFrom(settings);
    	m_proxypassword.loadSettingsFrom(settings);
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
    	m_gatk.validateSettings(settings);
    	m_use_dbsnp.validateSettings(settings);
    	m_dbsnp_file.validateSettings(settings);
    	m_use_interval.validateSettings(settings);
    	m_interval_file.validateSettings(settings);
    	m_variant_type.validateSettings(settings);
    	m_num_threads.validateSettings(settings);
    	
    	m_call_min_confidence.validateSettings(settings);
    	m_emit_min_confidence.validateSettings(settings);
    	m_pcr_error.validateSettings(settings);
    	m_contamination.validateSettings(settings);
    	m_heterozygosity.validateSettings(settings);
    	m_max_deletion_fraction.validateSettings(settings);
    	m_min_base_qual.validateSettings(settings);
    	m_indel_heterozygosity.validateSettings(settings);
    	m_min_indel_count.validateSettings(settings);
    	m_min_indel_frac.validateSettings(settings);
    	m_gap_open_pen.validateSettings(settings);
    	m_gap_cont_pen.validateSettings(settings);    	
    	m_baq.validateSettings(settings);
    	m_mbq.validateSettings(settings);
    	
    	//Proxy options
    	m_useproxy.validateSettings(settings);
    	m_proxyhost.validateSettings(settings);
    	m_proxyport.validateSettings(settings);
    	m_useproxyauth.validateSettings(settings);
    	m_proxyuser.validateSettings(settings);
    	m_proxypassword.validateSettings(settings);

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        
        // load internal data. 
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
       
        // save internal models. 
        // Everything written to output ports is saved automatically (data
        // returned by the execute method, models saved in the saveModelContent,
        // and user settings saved through saveSettingsTo - is all taken care 
        // of). Save here only the other internals that need to be preserved
        // (e.g. data used by the views).

    }

}


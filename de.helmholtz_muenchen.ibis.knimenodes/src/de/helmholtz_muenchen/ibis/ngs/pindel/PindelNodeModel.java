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
package de.helmholtz_muenchen.ibis.ngs.pindel;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

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
import org.knime.core.node.defaultnodesettings.SettingsModelDate;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.lofs.PathProcessor;


/**
 * This is the model implementation of Pindel.
 * 
 *
 * @author Marie-Sophie Friedl
 * @author Maximilian Hastreiter
 */
public class PindelNodeModel extends HTExecutorNodeModel {
    
    // the logger instance
    protected static final NodeLogger logger = NodeLogger.getLogger(PindelNodeModel.class);
        
    /** the settings key which is used to retrieve and 
        store the settings (from the dialog or from a settings file)    
       (package visibility to be usable from the dialog). */
    
    //path to pindel executable
    static final String CFGKEY_PINDEL="pindel";
    static final String DEF_PINDEL="";
	private final SettingsModelString m_pindel = new SettingsModelString(CFGKEY_PINDEL, DEF_PINDEL);
	
	//path to ref genome
	static final String CFGKEY_REFSEQFILE="refgenome";
	private final SettingsModelString m_refseqfile 			= new SettingsModelString(CFGKEY_REFSEQFILE,"");
	
	// restrict variant calling to a chromosome region
	static final String CFGKEY_INTERVAL="interval";
	static final boolean DEF_INTERVAL=false;
	private final SettingsModelBoolean m_interval = new SettingsModelBoolean(CFGKEY_INTERVAL, DEF_INTERVAL);
	
	// interval chromosome 
	static final String CFGKEY_CHROM ="chrom";
	static final String DEF_CHROM = "";
	private final SettingsModelString m_chrom = new SettingsModelString(CFGKEY_CHROM, DEF_CHROM);
	
	// interval start
	static final String CFGKEY_START="start";
	static final int DEF_START=0;
	static final int MIN_START=0;
	static final int MAX_START=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_start = new SettingsModelIntegerBounded(CFGKEY_START, DEF_START, MIN_START, MAX_START);
	
	// interval end
	static final String CFGKEY_END="end";
	static final int MIN_END=0;
	static final int MAX_END=Integer.MAX_VALUE;
	static final int DEF_END=250000000;
	private final SettingsModelInteger m_end = new SettingsModelIntegerBounded(CFGKEY_END, DEF_END, MIN_END, MAX_END);
	
	// path to config file
	static final String CFGKEY_CONFIG_FILE="config_file";
	static final String DEF_CONFIG_FILE="";
	private final SettingsModelString m_config_file= new SettingsModelString(CFGKEY_CONFIG_FILE, DEF_CONFIG_FILE);
	
	// create config file if previous node is Picard Tools CollectInsertSizeMetrics
	static final String CFGKEY_CREATE_CONFIG="create_config";
	static final boolean DEF_CREATE_CONFIG=false;
	private final SettingsModelBoolean m_create_config = new SettingsModelBoolean(CFGKEY_CREATE_CONFIG, DEF_CREATE_CONFIG);
	
	// convert pindel_D and pindel_SI to vcf
	static final String CFGKEY_VCF_OUT="vcf_out";
	static final boolean DEF_VCF_OUT=true;
	private final SettingsModelBoolean m_vcf_out = new SettingsModelBoolean(CFGKEY_VCF_OUT, DEF_VCF_OUT);
	
	// path to pindel2vcf converter
	static final String CFGKEY_VCF2PINDEL="vcf2pindel";
	static final String DEF_VCF2PINDEL="";
	private final SettingsModelString m_pindel2vcf = new SettingsModelString(CFGKEY_VCF2PINDEL, DEF_VCF2PINDEL);
	
	// number of threads
	static final String CFGKEY_THREADS="threads";
	static final int DEF_THREADS=1;
	static final int MIN_THREADS=1;
	static final int MAX_THREADS=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_threads = new SettingsModelIntegerBounded(CFGKEY_THREADS, DEF_THREADS, MIN_THREADS, MAX_THREADS);
	
	// bin size: divide reference in subsequences of X Mb
	static final String CFGKEY_BIN_SIZE="bin_size";
	static final int DEF_BIN_SIZE=10;
	static final int MIN_BIN_SIZE=1;
	static final int MAX_BIN_SIZE=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_bin_size = new SettingsModelIntegerBounded(CFGKEY_BIN_SIZE, DEF_BIN_SIZE, MIN_BIN_SIZE, MAX_BIN_SIZE);
	
	// minimum #mismatches for considering a read
	static final String CFGKEY_MIN_MATCH_BASES="min_match_bases";
	static final int DEF_MIN_MATCH_BASES=30;
	static final int MIN_MIN_MATCH_BASES=0;
	static final int MAX_MIN_MATCH_BASES=Integer.MAX_VALUE;
	private final SettingsModelInteger m_min_match_bases = new SettingsModelIntegerBounded(CFGKEY_MIN_MATCH_BASES, DEF_MIN_MATCH_BASES, MIN_MIN_MATCH_BASES, MAX_MIN_MATCH_BASES);
	
	// only alignment if no other position with less than this mismatches
	static final String CFGKEY_ADDITIONAL_MISMATCH="additional_mismatch";
	static final int DEF_ADDITIONAL_MISMATCH=1;
	static final int MIN_ADDITIONAL_MISMATCH=0;
	static final int MAX_ADDITIONAL_MISMATCH=Integer.MAX_VALUE;
	private final SettingsModelInteger m_additional_mismatch = new SettingsModelIntegerBounded(CFGKEY_ADDITIONAL_MISMATCH, DEF_ADDITIONAL_MISMATCH, MIN_ADDITIONAL_MISMATCH, MAX_ADDITIONAL_MISMATCH);
	
	// minimum matches requires at breakpoint
	static final String CFGKEY_MIN_MATCH_BP="min_match_breakpoint";
	static final int DEF_MIN_MATCH_BP=3;
	static final int MIN_MIN_MATCH_BP=0;
	static final int MAX_MIN_MATCH_BP=Integer.MAX_VALUE;
	private final SettingsModelInteger m_min_match_breakpoint = new SettingsModelIntegerBounded(CFGKEY_MIN_MATCH_BP, DEF_MIN_MATCH_BP, MIN_MIN_MATCH_BP, MAX_MIN_MATCH_BP);
	// sequencing error
	static final String CFGKEY_SEQ_ERR="seq_err";
	static final double DEF_SEQ_ERROR=0.05;
	static final double MIN_SEQ_ERROR=0;
	static final double MAX_SEQ_ERROR=1;
	private final SettingsModelDoubleBounded m_seq_err = new SettingsModelDoubleBounded(CFGKEY_SEQ_ERR, DEF_SEQ_ERROR, MIN_SEQ_ERROR, MAX_SEQ_ERROR);
	
	// only consider reads with less than this mismatches
	static final String CFGKEY_MAX_MISMATCH_RATE="max_mismatch_rate";
	static final double DEF_MAX_MISMATCH_RATE=0.1;
	static final double MIN_MAX_MISMATCH_RATE=0;
	static final double MAX_MAX_MISMATCH_RATE=1;
	private final SettingsModelDoubleBounded m_max_mismatch_rate = new SettingsModelDoubleBounded(CFGKEY_MAX_MISMATCH_RATE, DEF_MAX_MISMATCH_RATE, MIN_MAX_MISMATCH_RATE, MAX_MAX_MISMATCH_RATE);
	
	// reference name (filename)
	static final String CFGKEY_USE_REF_FILENAME="use_ref_filename";
	static final boolean DEF_USE_REF_FILENAME=true;
	private final SettingsModelBoolean m_use_ref_filename = new SettingsModelBoolean(CFGKEY_USE_REF_FILENAME, DEF_USE_REF_FILENAME);
	static final String CFGKEY_REFNAME="refname";
	static final String DEF_REFNAME="refname";
	private final SettingsModelString m_refname = new SettingsModelString(CFGKEY_REFNAME, DEF_REFNAME);
	
	// reference date (current date)
	static final String CFGKEY_USE_CUR_DATE="use_cur_date";
	static final boolean DEF_USE_CUR_DATE=true;
	private final SettingsModelBoolean m_use_cur_date = new SettingsModelBoolean(CFGKEY_USE_CUR_DATE, DEF_USE_CUR_DATE);
	static final String CFGKEY_REFDATE="refdate";
	private final SettingsModelDate m_refdate = new SettingsModelDate(CFGKEY_REFDATE);
	
	// minimum coverage reads (10)
	static final String CFGKEY_MIN_READS="min_reads";
	static final int DEF_MIN_READS=10;
	static final int MIN_MIN_READS=1;
	static final int MAX_MIN_READS=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_min_reads = new SettingsModelIntegerBounded(CFGKEY_MIN_READS, DEF_MIN_READS, MIN_MIN_READS, MAX_MIN_READS);
	
	// heterozygosity threshold (0.2)
	static final String CFGKEY_HETERO_FRAC="hetero_frac";
	static final double DEF_HETERO_FRAC=0.2;
	static final double MIN_HETERO_FRAC=0.0;
	static final double MAX_HETERO_FRAC=1.0;
	private final SettingsModelDoubleBounded m_hetero_frac = new SettingsModelDoubleBounded(CFGKEY_HETERO_FRAC, DEF_HETERO_FRAC, MIN_HETERO_FRAC, MAX_HETERO_FRAC);
	
	// homozygosity threshold (0.8)
	static final String CFGKEY_HOMO_FRAC="homo_frac";
	static final double DEF_HOMO_FRAC=0.8;
	static final double MIN_HOMO_FRAC=0.0;
	static final double MAX_HOMO_FRAC=1.0;
	private final SettingsModelDoubleBounded m_homo_frac = new SettingsModelDoubleBounded(CFGKEY_HOMO_FRAC, DEF_HOMO_FRAC, MIN_HOMO_FRAC, MAX_HOMO_FRAC);
	
	// GATK compatible
	static final String CFGKEY_GATK_COMP="gatk_comp";
	static final boolean DEF_GATK_COMP=true;
	private final SettingsModelBoolean m_gatk_comp = new SettingsModelBoolean(CFGKEY_GATK_COMP, DEF_GATK_COMP);
	
	// minimum reads supporting the variant (1)
	static final String CFGKEY_MIN_SUPP_READS="min_supp_reads";
	static final int DEF_MIN_SUPP_READS=1;
	static final int MIN_MIN_SUPP_READS=1;
	static final int MAX_MIN_SUPP_READS=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_min_supp_reads = new SettingsModelIntegerBounded(CFGKEY_MIN_SUPP_READS, DEF_MIN_SUPP_READS, MIN_MIN_SUPP_READS, MAX_MIN_SUPP_READS);
	
	// event supported on both strands? (false)
	static final String CFGKEY_BOTH_STRANDS="both_strands";
	static final boolean DEF_BOTH_STRANDS=false;
	private final SettingsModelBoolean m_both_strands = new SettingsModelBoolean(CFGKEY_BOTH_STRANDS, DEF_BOTH_STRANDS);
	
	// minimum size of variant (1)
	static final String CFGKEY_MIN_SIZE="min_size";
	static final int DEF_MIN_SIZE=1;
	static final int MIN_MIN_SIZE=1;
	static final int MAX_MIN_SIZE=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_min_size = new SettingsModelIntegerBounded(CFGKEY_MIN_SIZE, DEF_MIN_SIZE, MIN_MIN_SIZE, MAX_MIN_SIZE);
	
	// limit maximum size of variant
	static final String CFGKEY_LIMIT_SIZE="limit_size";
	static final boolean DEF_LIMIT_SIZE=false;
	private final SettingsModelBoolean m_limit_size = new SettingsModelBoolean(CFGKEY_LIMIT_SIZE, DEF_LIMIT_SIZE);
	
	// maximum size of variant (infinite)
	static final String CFGKEY_MAX_SIZE="max_size ";
	static final int DEF_MAX_SIZE=10;
	static final int MIN_MAX_SIZE=1;
	static final int MAX_MAX_SIZE=Integer.MAX_VALUE;
	private final SettingsModelIntegerBounded m_max_size = new SettingsModelIntegerBounded(CFGKEY_MAX_SIZE, DEF_MAX_SIZE, MIN_MAX_SIZE, MAX_MAX_SIZE);
	
	
	private int posBam;
//	private int posRef;
	
	// if previous node is CollectInsertSizeMetrics and position of ism file
	private boolean ISM=false;
	private int posISM;
	
	
    /**
     * Constructor for the node model.
     */
    protected PindelNodeModel() {
    
        // one incoming port and one outgoing port
    	super(1, 1);
        
        addSetting(m_pindel);
        addSetting(m_refseqfile);
        addSetting(m_interval);
        addSetting(m_chrom); 
        addSetting(m_start);
        addSetting(m_end);
        addSetting(m_config_file);
        addSetting(m_create_config);
        addSetting(m_vcf_out);
        addSetting(m_pindel2vcf);
        addSetting(m_threads);
        addSetting(m_bin_size);
        
        addSetting(m_min_match_bases);
        addSetting(m_additional_mismatch);
        addSetting(m_min_match_breakpoint);
        addSetting(m_seq_err);
        addSetting(m_max_mismatch_rate);
        
        addSetting(m_use_ref_filename);
        addSetting(m_refname);
        addSetting(m_use_cur_date);
        addSetting(m_refdate);
        addSetting(m_min_reads);
        addSetting(m_hetero_frac);
        addSetting(m_homo_frac);
        addSetting(m_gatk_comp);
        addSetting(m_min_supp_reads);
        addSetting(m_both_strands);
        addSetting(m_min_size);
        addSetting(m_limit_size);
        addSetting(m_max_size);
    	
        // disable dialog components from beginning on
    	m_chrom.setEnabled(false);
    	m_start.setEnabled(false);
    	m_end.setEnabled(false);
    	m_create_config.setEnabled(true);
    	
    	m_refname.setEnabled(false);
    	m_refdate.setEnabled(false);
    	m_max_size.setEnabled(false);
    	
    	}

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	
    	RunPindel rPindel = new RunPindel();
        //retrieve information from table
        DataRow r=inData[0].iterator().next();
                
        // check bam input file
        String inputfile=r.getCell(posBam).toString();
        
        // check if path is null
        if(inputfile.equals("")){
        	throw new Exception("No bam file available, something went wrong with the previous node!");
        }
        
        // check path to bam file
        if(!Files.exists(Paths.get(inputfile))){
        	throw new Exception("Path to input bam file: "+inputfile+" does not exist!");
        }
        
        // process path to input file
        String base=PathProcessor.getBase(inputfile);
        String fileextension=PathProcessor.getExt(inputfile);
        
        // check bam format
        if(!fileextension.equals("bam")){
        	throw new Exception("Input file is not in bam format!");
        }
        
        // check bam file index
        if(!Files.exists(Paths.get(inputfile+".bai"))){
        	
	        if(!Files.exists(Paths.get(base+".bai"))){
	        		throw new Exception("Missing bam file index: "+base+".bai");
	        }
	        else{
	        	logger.info("Renaming index file");
	        	Files.copy(Paths.get(base+".bai"), Paths.get(inputfile+".bai"));
	        }
        }
        
        // reference file
        
        String reffile=m_refseqfile.getStringValue();
        
        // path to reffile should not be null
        if(reffile.equals("")){
        	throw new Exception("No reference file available!");
        }
        
        // check path to reference path
        if(!Files.exists(Paths.get(reffile))){
        	throw new Exception("Reference sequence file: "+reffile+" does not exist!");
        }
        
        // ism metrics data
        String ismfile="";
        if(ISM && m_create_config.getBooleanValue()){
        	ismfile=r.getCell(posISM).toString();
        	
        	if(ismfile.equals("")){
        		throw new Exception("Something went wrong: Missing Insert Size Metrics File from previous node");
        	}
        	
        	if(!Files.exists(Paths.get(ismfile))){
        		throw new Exception("Insert Size Metrics File: "+ismfile+" does not exist!");
        	}
        }
        
        // check node options
        
        // check path to pindel file 
        if(m_pindel.getStringValue().equals("")){
        	throw new Exception("Missing path to pindel executable: You have to configure the node before executing it!");
        }
        
        String pindelfile = m_pindel.getStringValue();
        
        if(!Files.exists(Paths.get(pindelfile))){
        	throw new Exception("Path to pindel executable: "+pindelfile+" does not exist");
        }
        
       
        //check interval
        String interval="";
        if(m_interval.getBooleanValue()){
        	interval=m_chrom.getStringValue()+":"+m_start.getIntValue()+"-"+m_end.getIntValue();
        }

        
        //check pindel2vcf
        if(m_vcf_out.getBooleanValue()){
        	if(m_pindel2vcf.equals("")){
        		throw new Exception("Missing path to pindel2vcf converter: You have to configure the node before executing it!");
        	}
        	if(!Files.exists(Paths.get(m_pindel2vcf.getStringValue()))){
        		throw new Exception("Path to pindel2vcf converter: "+m_pindel2vcf.getStringValue()+" does not exist");
        	}
        }
        
        String configfile= "";	
        //create config file
        if(ISM && m_create_config.getBooleanValue()){
        	
        	configfile=PathProcessor.createOutputFile(base, "config", "pindel");
        	rPindel.PindelConfig(inputfile, ismfile, configfile);

        }else if(!m_create_config.getBooleanValue()){
        	configfile= m_config_file.getStringValue();
        }

        
        //run Pindel
        String pout= rPindel.createOutputFilePindel(base, "pindel", "D");
        
        int [] resources = new int[] {m_threads.getIntValue(), m_bin_size.getIntValue()};
        double [] sen_spec = new double[] {m_min_match_bases.getIntValue(), m_additional_mismatch.getIntValue(), m_min_match_breakpoint.getIntValue(), m_seq_err.getDoubleValue(), m_max_mismatch_rate.getDoubleValue()};

        
        rPindel.Pindel(exec, pindelfile, configfile, reffile, pout, interval, resources, sen_spec);
        
        //check output files: pindel_D deletions, pindel_SI small insertions
        String pindeldeletions=pout+"_D";
        if(!Files.exists(Paths.get(pindeldeletions))){
        	throw new Exception("Something went wrong executing pindel: misssing output file "+pindeldeletions);
        }
        String pindelinsertions=pout+"_SI";
        if(!Files.exists(Paths.get(pindelinsertions))){
        	throw new Exception("Something went wrong executing pindel: misssing output file "+pindelinsertions);
        }

        
        //convert pindel output to vcf
    	String delout=pindeldeletions+".vcf";
    	String inout=pindelinsertions+".vcf";

    	
        if(m_vcf_out.getBooleanValue()){
        	
        	String refname="";
        	if(m_use_ref_filename.getBooleanValue()){
        		refname=Paths.get(PathProcessor.getBase(reffile)).getFileName().toString();
        	}
        	else{
        		refname=m_refname.getStringValue();
        	}
        	
        	DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd");
        	
        	Date date=null;
        	if(m_use_cur_date.getBooleanValue()){
        		date = new Date();
        	}
        	else{
        		date = m_refdate.getDate();
        	}
        	String refdate = dateFormat.format(date);
        	
        	boolean [] flags = new boolean []{m_gatk_comp.getBooleanValue(), m_both_strands.getBooleanValue(), m_limit_size.getBooleanValue()};
        	double [] numbers = new double []{m_min_reads.getIntValue(), m_hetero_frac.getDoubleValue(), m_homo_frac.getDoubleValue(), m_min_supp_reads.getIntValue(), m_min_size.getIntValue(), m_max_size.getIntValue()};
 
        	rPindel.Pindel2VCF(exec, m_pindel2vcf.getStringValue(), reffile, refname, refdate, pindeldeletions, delout, flags, numbers);
        	rPindel.Pindel2VCF(exec, m_pindel2vcf.getStringValue(), reffile, refname, refdate, pindelinsertions, inout, flags, numbers);
    	}
        
        //create output table
        
    	// create column specifications
    	DataColumnSpec [] colspec = new DataColumnSpec[2];
    	if(m_vcf_out.getBooleanValue()){
	    	colspec[0]=new DataColumnSpecCreator("Path2VCFdeletionsFile", StringCell.TYPE).createSpec();
	    	colspec[1]=new DataColumnSpecCreator("Path2VCFinsertionsFile", StringCell.TYPE).createSpec();    		
    	}
    	else{
	    	colspec[0]=new DataColumnSpecCreator("Path2PindelDFile", StringCell.TYPE).createSpec();
	    	colspec[1]=new DataColumnSpecCreator("Path2PindelSIFile", StringCell.TYPE).createSpec();
    	}
  	
    	//create table
	    DataTableSpec outspec=new DataTableSpec(colspec);
	    BufferedDataContainer c = exec.createDataContainer(outspec);
	    
	    // fill string cells
	    DataCell [] row = new DataCell [2];
	    if(m_vcf_out.getBooleanValue()){
	    	row[0]=new StringCell(delout);
	    	row[1]=new StringCell(inout);	    	
	    }
	    else{
	    	row[0]=new StringCell(pindeldeletions);
    		row[1]=new StringCell(pindelinsertions);
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
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	//check interval
    	if (m_interval.getBooleanValue()){
    		//check chromosome name
    		if(m_chrom.getStringValue().equals("")){
    			throw new InvalidSettingsException("You have to enter a chromosome name");
    		}
    		//check if start>=end
        	if(m_end.getIntValue()<m_start.getIntValue()){
        		throw new InvalidSettingsException("Invalid interval "+m_chrom.getStringValue()+":"+m_start.getIntValue()+"-"+m_end.getIntValue()+": end < start!");
        	}
    	}
    	
    	if(m_hetero_frac.getDoubleValue()>= m_homo_frac.getDoubleValue()){
    		throw new InvalidSettingsException("Fraction for heterozygosity has to be smaller than fraction for homozygosity!");
    	}
    	
    	if(m_min_size.getIntValue()> m_max_size.getIntValue()){
    		throw new InvalidSettingsException("Minimum variant size has to be samller than maximum variant size");
    	}
    	    	
    	String [] cols = inSpecs[0].getColumnNames();
    	    	
    	//check if previous node is CollectInsertSizeMetrics
    	if(inSpecs[0].containsName("Path2ISMetrics")){
    		logger.info("Previous Node is PicardTools: CollectInsertSizeMetrics");
    		ISM=true;
    		m_create_config.setEnabled(true);
    		
        	//get positions of metrics file of in port table
        	for (int i=0; i<cols.length; i++){
        		if(cols[i].equals("Path2ISMetrics")){
        			posISM=i;
        		}
        	}
    	}
    	else
    	{
    		ISM=false;
    		m_create_config.setEnabled(false);
    		m_config_file.setEnabled(true);
    	}
    	
        //check pindel config file
        String configfile="";
        
        if(!ISM && m_create_config.getBooleanValue()){
        	throw new InvalidSettingsException("Previous node is not Picard Tools: CollectInsertSizeMetrics! You have to run in order to create a config file.");
        }
        else if (!ISM || !m_create_config.getBooleanValue()){
	        if(m_config_file.getStringValue().equals("")){
	        	throw new InvalidSettingsException("Missing path to pindel config file: You have to configure the node before executing it!");
	        }

	        configfile= m_config_file.getStringValue();
	        
	        if(!Files.exists(Paths.get(configfile))){
	        	throw new InvalidSettingsException("Path to pindel config: "+configfile+" does not exist");
	        }
        }
    	
    	//Check if bam file is available
    	posBam=CompatibilityChecker.getFirstIndexCellType(inSpecs[0], "BAMCell");
    	if(posBam==-1){
    		throw new InvalidSettingsException("No BAM found in input table. Wrong node/input connected?");
    	}
    	    	
    	// create column specifications
    	DataColumnSpec [] colspec = new DataColumnSpec[2];
    	if(m_vcf_out.getBooleanValue()){
	    	colspec[0]=new DataColumnSpecCreator("Path2VCFdeletionsFile", VCFCell.TYPE).createSpec();
	    	colspec[1]=new DataColumnSpecCreator("Path2VCFinsertionsFile", VCFCell.TYPE).createSpec();    		
    	}
    	else{
	    	colspec[0]=new DataColumnSpecCreator("Path2PindelDFile", StringCell.TYPE).createSpec();
	    	colspec[1]=new DataColumnSpecCreator("Path2PindelSIFile", StringCell.TYPE).createSpec();
    	}
  	
    	//create table
	    DataTableSpec outspec=new DataTableSpec(colspec);
    	
    	
        return new DataTableSpec[]{outspec};
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


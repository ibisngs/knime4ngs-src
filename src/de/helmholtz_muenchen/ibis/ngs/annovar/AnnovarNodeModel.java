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

package de.helmholtz_muenchen.ibis.ngs.annovar;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of Annovar.
 * 
 *
 * @author Sebastian Kopetzky
 */
public class AnnovarNodeModel extends NodeModel {
	/**
	 *             Arguments to download databases or perform annotations
                --downdb                    download UCSC Genome Browser annotation database
                --geneanno                  annotate variants by functional consequences on genes
                --regionanno                annotate variants by targetting specific genomics regions
                --filter                    filter variants based on a position list
                --webfrom <string>          specify the source of database (default usually works fine)
        
            Arguments to control input and output
                --outfile <file>            output file prefix
                --zerostart                 input query file uses half-open zero-start coordinate
                --dbtype <string>           database type
                --buildver <string>         genome build version (default: hg18 for human)
                --gff3dbfile <file>         specify the GFF3 DB file used in region-based annotation
                --genericdbfile <file>      specify the generic DB file used in filter-based annotation
                --vcfdbfile <file>          specify the DB file in VCF format in filter-based annotation
                --bedfile <file>            specify a BED file in region-based annotation
                --time                      print out local time during program run
                --separate                  separately print out all function of a variant (default: one line per variant)
                --colsWanted <string>       specify which columns to output in -regionanno by comma-delimited numbers
                --comment                   print out comment line (those starting with #) in output files 
                --scorecolumn <int>         the column with scores in database file (for region-based annotation)
                --exonsort                  sort the exon number in output line (for gene-based annotation)
                --transcript_function       use transcript name rather than gene name in gene-based annotation output
                --hgvs                      use HGVS format for exonic annotation (c.122C>T rather than c.C122T)
                --otherinfo                 in filter-based annotation, print out additional columns in database file
                --infoasscore               in filter-based annotation, use INFO field in VCF file as score in output
                --seq_padding               if set, create a new file with cDNA sequence padded by this much either side
        
            Arguments to fine-tune the annotation procedure
                --batchsize <int>           batch size for processing variants per batch (default: 5m)
                --genomebinsize <int>       bin size to speed up search (default: 100k for -geneanno, 10k for -regionanno)
                --expandbin <int>           check nearby bin to find neighboring genes (default: 2m/genomebinsize)
                --neargene <int>            distance threshold to define upstream/downstream of a gene
                --score_threshold <float>   minimum score of DB regions to use in annotation
                --reverse                   reverse directionality to compare to score_threshold
                --normscore_threshold <float> minimum normalized score of DB regions to use in annotation
                --rawscore                  output includes the raw score (not normalized score) in UCSC Browser Track
                --minqueryfrac <float>      minimum percentage of query overlap to define match to DB (default: 0)
                --splicing_threshold <int>  distance between splicing variants and exon/intron boundary (default: 2)
                --indel_splicing_threshold <int>    if set, use this value for allowed indel size for splicing variants (default: --splicing_threshold)
                --maf_threshold <float>     filter 1000G variants with MAF above this threshold (default: 0)
                --sift_threshold <float>    SIFT threshold for deleterious prediction (default: 0.05)
                --precedence <string>       comma-delimited to specify precedence of variant function (default: exonic>intronic...)
                --indexfilter_threshold <float>     controls whether filter-based annotation use index if this fraction of bins need to be scanned (default: 0.9)
       
           Arguments to control memory usage
                --memfree <int>             ensure minimum amount of free system memory (default: 100000, in the order of kb)
                --memtotal <int>            limit total amount of memory used by ANNOVAR (default: 0, unlimited, in the order of kb)
                --chromosome <string>       examine these specific chromosomes in database file

	 */
	
	private final NodeLogger LOGGER = getLogger();
	
	
	/**
	 * Input arguments
	 */
	public static final String CFGKEY_QUERYFILE="queryfile";
	public static final String CFGKEY_TABLENAME="tablename";
	public static final String CFGKEY_DATABASELOCATION="databaselocation";
	public static final String CFGKEY_PATH2ANNOVAR="path2annovar";
	
	/**
	 * Input checkboxes
	 */
	public static final String CFGKEY_USEQUERYFILE="usequeryfile";
	public static final String CFGKEY_USETABLENAME="usetablename";
	
	/**
	 *  Arguments to download databases or perform annotations
	 */
	public static final String CFGKEY_METHOD="method"; //geneanno , regionanno, filter
	public static final String CFGKEY_WEBFROM="webfrom";
	/**
	 *  Arguments to control input and output
	 */
	public static final String CFGKEY_OUTFILE="outfile";
	public static final String CFGKEY_DBTYPE="dbtype";
	public static final String CFGKEY_BUILDVER="buildver";
	public static final String CFGKEY_GFF3DBFILE="gff3dbfile";
	public static final String CFGKEY_GENERICDBFILE="genericdbfile";
	public static final String CFGKEY_VCFDBFILE="vcfdbfile";
	public static final String CFGKEY_BEDFILE="bedfile";
	public static final String CFGKEY_SEPARATE="separate";
	public static final String CFGKEY_COLSWANTED="colswanted";
	public static final String CFGKEY_COMMENT="comment";
	public static final String CFGKEY_SCORECOLUMN="scorecolumn";
	public static final String CFGKEY_EXONSORT="exonsort";
	public static final String CFGKEY_TRANSCRIPTFUNCTION="transcriptfunction";
	public static final String CFGKEY_HGVS="hgvs";
	public static final String CFGKEY_OTHERINFO="otherinfo";
	public static final String CFGKEY_INFOASSCORE="infoasscore";
	public static final String CFGKEY_SEQPADDING="seqpadding";
	
	/**
	 * In/Out Checkboxes
	 */
//	public static final String CFGKEY_USEOUTFILE="useoutfile";
//	public static final String CFGKEY_USEDBTYPE="usedbtype";
//	public static final String CFGKEY_USEBUILDVER="usebuildver";
	public static final String CFGKEY_USEGFF3DBFILE="usegff3dbfile";
	public static final String CFGKEY_USEGENERICDBFILE="usegenericdbfile";
	public static final String CFGKEY_USEVCFDBFILE="usevcfdbfile";
	public static final String CFGKEY_USEBEDFILE="usebedfile";
	public static final String CFGKEY_USECOLSWANTED="usecolsWanted";
	public static final String CFGKEY_USESCORECOLUMN="usescorecolumn";

	/**
	 * Arguments to fine-tune the annotation procedure
	 */
	public static final String CFGKEY_BATCHSIZE="batchsize";
	public static final String CFGKEY_GENOMEBINSIZE="genomebinsize";
	public static final String CFGKEY_EXPANDBIN="expandbin";
	public static final String CFGKEY_NEARGENE="neargene";
	public static final String CFGKEY_SCORETHRESHOLD="scorethreshold";
	public static final String CFGKEY_REVERSE="reverse";
	public static final String CFGKEY_NORMSCORETHRESHOLD="normscorethreshold";
	public static final String CFGKEY_RAWSCORE="rawscore";
	public static final String CFGKEY_MINQUERYFRAC="minqueryfrac";
	public static final String CFGKEY_SPLICINGTHRESHOLD="splcingthreshold";
	public static final String CFGKEY_INDELSPLICINGTHRESHOLD="indelsplicingthreshold";
	public static final String CFGKEY_MAFTHRESHOLD="mafthreshold";
	public static final String CFGKEY_SIFTTHRESHOLD="siftthreshold";
	public static final String CFGKEY_PRECEDENCE="precedence";
	public static final String CFGKEY_INDEXFILTERTHRESHOLD="indexthreshold";
	
	/**
	 * Fine-tune checkboxes
	 */
	public static final String CFGKEY_USENEARGENE="useneargene";
	public static final String CFGKEY_USESCORETHRESHOLD="usescorethreshold";
	public static final String CFGKEY_USENORMSCORETHRESHOLD="usenormscorethreshold";
	public static final String CFGKEY_USESPLICINGTHRESHOLD="usesplicingthreshold";
	public static final String CFGKEY_USEINDELSPLICINGTHRESHOLD="useindelsplicingthreshold";
	public static final String CFGKEY_USESIFTTHRESHOLD="usesiftthreshold";
	
	
	/**
	 * Arguments to control memory usage
	 */
	public static final String CFGKEY_MEMFREE="memfree";
	public static final String CFGKEY_MEMTOTAL="memtotal";
	public static final String CFGKEY_CHROMOSOME="chromosome";
	
	/**
	 * Memory checkbox
	 */
	public static final String CFGKEY_USECHROMOSOME="usechromosome";
	
	/**
	 * Arguments to fine-tune the annotation procedure Defaults
	 */
	public static final int DEFAULT_BATCHSIZE=5;
	public static final int DEFAULT_GENOMEBINSIZE=100;
	public static final int DEFAULT_EXPANDBIN=20;
	public static final double DEFAULT_MINQUERYFRAC=0;
	public static final int DEFAULT_SPLICINGTHRESHOLD=2;
	public static final double DEFAULT_MAFTHRESHOLD=0;
	public static final double DEFAULT_SIFTTHRESHOLD=0.05;
	public static final String DEFAULT_PRECEDENCE="exonic,intronic";
	public static final double DEFAULT_INDEXFILTERTHRESHOLD=0.9;
	/**
	 * Arguments to control memory usage Defaults
	 */
	public static final int DEFAULT_MEMFREE=100000;
	public static final int DEFAULT_MEMTOTAL=0;
	
	
	/**
	 * Input arguments
	 */
	private final SettingsModelString m_queryfile = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_QUERYFILE,"");
	private final SettingsModelString m_tablename = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_TABLENAME,"");
	private final SettingsModelString m_databaselocation = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_DATABASELOCATION,"");
	private final SettingsModelString m_path2annovar = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_PATH2ANNOVAR,"");
	/**
	 * Input checkboxes
	 */
//	private final SettingsModelBoolean m_usequeryfile = new SettingsModelBoolean(
//			AnnovarNodeModel.CFGKEY_USEQUERYFILE,false);
	private final SettingsModelBoolean m_usetablename = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USETABLENAME,false);
	
	/**
	 * Arguments to download databases or perform annotations
	 */
	private final SettingsModelString m_method = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_METHOD,"");
//	private final SettingsModelString m_webfrom = new SettingsModelString(
//			AnnovarNodeModel.CFGKEY_WEBFROM,"");
	/**
	 *  Arguments to control input and output models
	 */
	private final SettingsModelString m_outfile = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_OUTFILE,"");
	private final SettingsModelString m_dbtype = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_DBTYPE,"");
	private final SettingsModelString m_buildver = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_BUILDVER,"");
	private final SettingsModelString m_gff3dbfile = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_GFF3DBFILE,"");
	private final SettingsModelString m_genericdbfile = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_GENERICDBFILE,"");
	private final SettingsModelString m_vcfdbfile = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_VCFDBFILE,"");
	private final SettingsModelString m_bedfile = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_BEDFILE,"");
	private final SettingsModelBoolean m_separate = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_SEPARATE, false);
	private final SettingsModelString m_colswanted = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_COLSWANTED,"");
	private final SettingsModelBoolean m_comment = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_COMMENT, false);
	private final SettingsModelIntegerBounded m_scorecolumn = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_SCORECOLUMN,0,0,Integer.MAX_VALUE);
	private final SettingsModelBoolean m_exonsort = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_EXONSORT, false);
	private final SettingsModelBoolean m_transcriptfunction = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_TRANSCRIPTFUNCTION, false);
	private final SettingsModelBoolean m_hgvs = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_HGVS, false);
	private final SettingsModelBoolean m_otherinfo = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_OTHERINFO, false);
	private final SettingsModelBoolean m_infoasscore = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_INFOASSCORE, false);
	private final SettingsModelBoolean m_seqpadding = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_SEQPADDING, false);
	
	/**
	 * In/Out Checkboxes
	 */
//	private final SettingsModelBoolean m_useoutfile = new SettingsModelBoolean(
//			AnnovarNodeModel.CFGKEY_USEOUTFILE, false);
//	private final SettingsModelBoolean m_usedbtype = new SettingsModelBoolean(
//			AnnovarNodeModel.CFGKEY_USEDBTYPE, false);
//	private final SettingsModelBoolean m_usebuildver = new SettingsModelBoolean(
//			AnnovarNodeModel.CFGKEY_USEBUILDVER, false);
	private final SettingsModelBoolean m_usegff3dbfile = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USEGFF3DBFILE, false);
	private final SettingsModelBoolean m_usegenericdbfile = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USEGENERICDBFILE, false);
	private final SettingsModelBoolean m_usevcfdbfile = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USEVCFDBFILE, false);
	private final SettingsModelBoolean m_usebedfile = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USEBEDFILE, false);
	private final SettingsModelBoolean m_usecolswanted = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USECOLSWANTED, false);
	private final SettingsModelBoolean m_usescorecolumn = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USESCORECOLUMN, false);
	

	/**
	 * Arguments to fine-tune the annotation procedure
	 */
	private final SettingsModelIntegerBounded m_batchsize = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_BATCHSIZE,DEFAULT_BATCHSIZE,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_genomebinsize = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_GENOMEBINSIZE,DEFAULT_GENOMEBINSIZE,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_expandbin = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_EXPANDBIN,DEFAULT_EXPANDBIN,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_neargene = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_NEARGENE,0,0,Integer.MAX_VALUE);
	private final SettingsModelDoubleBounded m_scorethreshold = new SettingsModelDoubleBounded(
			AnnovarNodeModel.CFGKEY_SCORETHRESHOLD,0,0,Double.MAX_VALUE);
	private final SettingsModelBoolean m_reverse = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_REVERSE, false);
	private final SettingsModelIntegerBounded m_normscorethreshold = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_NORMSCORETHRESHOLD,0,0,Integer.MAX_VALUE);
	private final SettingsModelBoolean m_rawscore = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_RAWSCORE, false);
	private final SettingsModelDoubleBounded m_minqueryfrac = new SettingsModelDoubleBounded(
			AnnovarNodeModel.CFGKEY_MINQUERYFRAC,DEFAULT_MINQUERYFRAC,0,Double.MAX_VALUE);
	private final SettingsModelIntegerBounded m_splicingthreshold = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_SPLICINGTHRESHOLD,DEFAULT_SPLICINGTHRESHOLD,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_indelsplicingthreshold = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_INDELSPLICINGTHRESHOLD,0,0,Integer.MAX_VALUE);
	private final SettingsModelDoubleBounded m_mafthreshold = new SettingsModelDoubleBounded(
			AnnovarNodeModel.CFGKEY_MAFTHRESHOLD,DEFAULT_MAFTHRESHOLD,0,Double.MAX_VALUE);
	private final SettingsModelDoubleBounded m_siftthreshold = new SettingsModelDoubleBounded(
			AnnovarNodeModel.CFGKEY_SIFTTHRESHOLD,DEFAULT_SIFTTHRESHOLD,0,Double.MAX_VALUE);
	private final SettingsModelString m_precedence = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_PRECEDENCE,DEFAULT_PRECEDENCE);
	private final SettingsModelDoubleBounded m_indexfilterthreshold = new SettingsModelDoubleBounded(
			AnnovarNodeModel.CFGKEY_INDEXFILTERTHRESHOLD,DEFAULT_INDEXFILTERTHRESHOLD,0,Double.MAX_VALUE);
	
	/**
	 * Fine-tune checkboxes
	 */
	private final SettingsModelBoolean m_useneargene = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USENEARGENE, false);
	private final SettingsModelBoolean m_usescorethreshold = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USESCORETHRESHOLD, false);
	private final SettingsModelBoolean m_usenormscorethreshold = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USENORMSCORETHRESHOLD, false);
	private final SettingsModelBoolean m_usesplicingthreshold = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USESPLICINGTHRESHOLD, false);
	private final SettingsModelBoolean m_useindelsplicingthreshold = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USEINDELSPLICINGTHRESHOLD, false);
	private final SettingsModelBoolean m_usesiftthreshold = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USESIFTTHRESHOLD, false);
	
	
	/**
	 * Arguments to control memory usage models
	 */
	private final SettingsModelIntegerBounded m_memfree = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_MEMFREE,DEFAULT_MEMFREE,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_memtotal = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_MEMTOTAL,DEFAULT_MEMTOTAL,0,Integer.MAX_VALUE);
	private final SettingsModelString m_chromosome = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_CHROMOSOME,"");
	/**
	 * Memory checkbox
	 */
	private final SettingsModelBoolean m_usechromosome = new SettingsModelBoolean(
			AnnovarNodeModel.CFGKEY_USECHROMOSOME, false);
	

	private boolean optionalPort=false;
	
	
    /**
     * Constructor for the node model.
     */
    protected AnnovarNodeModel() {
    
       super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(0));

        
       m_neargene.setEnabled(false);
       m_scorethreshold.setEnabled(false);
       m_normscorethreshold.setEnabled(false);
       m_splicingthreshold.setEnabled(false);
       m_indelsplicingthreshold.setEnabled(false);     
       m_chromosome.setEnabled(false);
       //m_outfile.setEnabled(false);
      // m_dbtype.setEnabled(false);
       m_gff3dbfile.setEnabled(false);
       m_genericdbfile.setEnabled(false);
       m_vcfdbfile.setEnabled(false);
       m_bedfile.setEnabled(false);
       m_colswanted.setEnabled(false);
       m_scorecolumn.setEnabled(false);
       m_otherinfo.setEnabled(false);
       m_infoasscore.setEnabled(false);
       m_tablename.setEnabled(false);     
       m_mafthreshold.setEnabled(false);
       m_siftthreshold.setEnabled(false);
       m_indexfilterthreshold.setEnabled(false);
       m_minqueryfrac.setEnabled(false);
       
       m_usegff3dbfile.setEnabled(false);
       m_usevcfdbfile.setEnabled(false);
       m_usegenericdbfile.setEnabled(false);
       m_usebedfile.setEnabled(false);
       m_usescorecolumn.setEnabled(false);
       m_usecolswanted.setEnabled(false);
       
       m_normscorethreshold.setEnabled(false);
       m_scorethreshold.setEnabled(false);
       m_usenormscorethreshold.setEnabled(false);
       m_usescorethreshold.setEnabled(false);
       m_usesiftthreshold.setEnabled(false);
       m_reverse.setEnabled(false);
       
      
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	/**Initialize logfile**/
    	String folder = m_outfile.getStringValue()+"/";
    	String logfile = folder +"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("Annovar"));
    	/**logfile initialized**/
       	//annotate_variation.pl [arguments] <query-file|table-name> <database-location>
    	ArrayList<String> command = new ArrayList<String>();
    	/**
    	 * Input arguments
    	 */

    	if(optionalPort){
    		command.add(inData[0].iterator().next().getCell(1).toString()+"/annotate_variation.pl");
    	}
    	else{
    		command.add(m_path2annovar.getStringValue()+"/annotate_variation.pl");
    	}
    	
    	
    	if(m_method.getStringValue().equals("geneanno")){command.add("--geneanno");}
    	if(m_method.getStringValue().equals("regionanno")){command.add("--regionanno");}
    	if(m_method.getStringValue().equals("filter")){command.add("--filter");}
    	
    	/**
    	 *  Arguments to control input and output models
    	 */
    	String outfile=m_outfile.getStringValue();
    	if(m_queryfile.isEnabled()){
    		String qfile = m_queryfile.getStringValue();
    		qfile=qfile.substring(qfile.lastIndexOf("/"), qfile.lastIndexOf("."));
    		outfile+=qfile;
    	}else{
    			String tname=m_tablename.getStringValue();
    			outfile+="/"+tname;
    	}

    	command.add("--outfile "+outfile);
    	command.add("--dbtype "+m_dbtype.getStringValue());
    	command.add("--buildver "+m_buildver.getStringValue());
    	if(m_usegff3dbfile.getBooleanValue()){command.add("--gff3dbfile "+m_gff3dbfile.getStringValue());}
    	if(m_usegenericdbfile.getBooleanValue()){command.add("--genericdbfile "+m_genericdbfile.getStringValue());}
    	if(m_usevcfdbfile.getBooleanValue()){command.add("--vcfdbfile "+m_vcfdbfile.getStringValue());}
    	if(m_usebedfile.getBooleanValue()){command.add("--bedfile "+m_bedfile.getStringValue());}
    	if(m_separate.getBooleanValue()){command.add("--separate");}
    	if(m_usecolswanted.getBooleanValue()){command.add("--colsWanted "+m_colswanted.getStringValue());}
    	if(m_comment.getBooleanValue()){command.add("--comment");}
    	if(m_usescorecolumn.getBooleanValue()){command.add("--scorecolumn "+m_scorecolumn.getIntValue());}
    	if(m_exonsort.getBooleanValue()){command.add("--exonsort");}
    	if(m_transcriptfunction.getBooleanValue()){command.add("--transcript_function");}
    	if(m_hgvs.getBooleanValue()){command.add("--hgvs");}
      	if(m_otherinfo.getBooleanValue()){command.add("--otherinfo");}
      	if(m_infoasscore.getBooleanValue()){command.add("--infoasscore");}
      	if(m_seqpadding.getBooleanValue()){command.add("--seq_padding");}
      	
    	/**
    	 * Arguments to fine-tune the annotation procedure
    	 */
      	command.add("--batchsize "+m_batchsize.getIntValue()+"m");
      	command.add("--genomebinsize "+m_genomebinsize.getIntValue()+"k");
      	command.add("--expandbin "+m_expandbin.getIntValue());
      	if(m_useneargene.getBooleanValue()){command.add("--neargene "+m_neargene.getIntValue());}
      	if(m_usescorethreshold.getBooleanValue()){command.add("--score_threshold "+m_scorethreshold.getDoubleValue());}
      	if(m_reverse.getBooleanValue()){command.add("--reverse");}
      	if(m_usenormscorethreshold.getBooleanValue()){command.add("--normscore_threshold "+m_normscorethreshold.getIntValue());}
      	if(m_rawscore.getBooleanValue()){command.add("--rawscore");}
      	if(m_minqueryfrac.isEnabled()){command.add("--minqueryfrac "+m_minqueryfrac.getDoubleValue());}
      	if(m_usesplicingthreshold.getBooleanValue()){command.add("--splicing_threshold "+m_splicingthreshold.getIntValue());}
      	if(m_useindelsplicingthreshold.getBooleanValue()){command.add("--indel_splicing_threshold "+m_indelsplicingthreshold.getIntValue());}
      	if(m_mafthreshold.isEnabled()){command.add("--maf_threshold "+m_mafthreshold.getDoubleValue());}
      	if(m_siftthreshold.isEnabled()){command.add("--sift_threshold "+m_siftthreshold.getDoubleValue());}
      	if(m_precedence.isEnabled()){command.add("--precedence "+m_precedence.getStringValue());}
      	if(m_indexfilterthreshold.isEnabled()){command.add("--indexfilter_threshold "+m_indexfilterthreshold.getDoubleValue());}

    	/**
    	 * Arguments to control memory usage models
    	 */
      	command.add("--memfree "+m_memfree.getIntValue());
      	command.add("--memtotal "+m_memtotal.getIntValue());
      	if(m_chromosome.isEnabled()){command.add("--chromosome "+m_chromosome.getStringValue());}
      	
      	/**
      	 * Output files
      	 */
		if(optionalPort){ //get name from inData array
			String qfile = inData[0].iterator().next().getCell(0).toString();
			command.add(qfile);
			command.add(m_databaselocation.getStringValue()); 
		}else{
			if(m_usetablename.getBooleanValue())
			{
				command.add(m_tablename.getStringValue());
			}else{
				command.add(m_queryfile.getStringValue());
				command.add(m_databaselocation.getStringValue()); 
			}
		}

    	/**
    	 * Execute
    	 */
    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
    	logBuffer.append(ShowOutput.getNodeEndTime());
    	ShowOutput.writeLogFile(logBuffer);
    	
    	
        return new BufferedDataTable[]{};
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

    	
    	
    	
		//Check OptionalInputPort
		try{
			inSpecs[0].getColumnNames();
			optionalPort=true;
			String[] colnames = inSpecs[0].getColumnNames();
			if(colnames[0].equals("Path2OutFile") && colnames[1].equals("InstallPath")){
				m_queryfile.setEnabled(false);
	    		m_usetablename.setEnabled(false);
	    		m_path2annovar.setEnabled(false);
			}
			else{
				throw new InvalidSettingsException("The in-port is invalid! Don't use an in-port or connect the right node.");
			}

			//m_bfast.setEnabled(false);
    	}catch(NullPointerException npe){
    		m_queryfile.setEnabled(true);
    		m_usetablename.setEnabled(true);
    		m_path2annovar.setEnabled(true);
    	}
    	/*if(inSpecs[0] != null){
    		OptionalPort=true;
    	}
    	else{
    		m_queryfile.setEnabled(false);
    		m_usetablename.setEnabled(false);
    	}*/
    	
    	
    	if(m_siftthreshold.isEnabled() && !m_dbtype.getStringValue().equals("avsift")){
            throw new InvalidSettingsException("SIFT threshold can only be used for dbtype=avsift!");
    	}
        return new DataTableSpec[]{};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_batchsize.saveSettingsTo(settings);
         m_bedfile.saveSettingsTo(settings);
         m_buildver.saveSettingsTo(settings);
         m_chromosome.saveSettingsTo(settings);
         m_colswanted.saveSettingsTo(settings);
         m_comment.saveSettingsTo(settings);
         m_databaselocation.saveSettingsTo(settings);
         m_dbtype.saveSettingsTo(settings);
         m_exonsort.saveSettingsTo(settings);
         m_expandbin.saveSettingsTo(settings);
         m_genericdbfile.saveSettingsTo(settings);
         m_genomebinsize.saveSettingsTo(settings);
         m_gff3dbfile.saveSettingsTo(settings);
         m_hgvs.saveSettingsTo(settings);
         m_indelsplicingthreshold.saveSettingsTo(settings);
         m_indexfilterthreshold.saveSettingsTo(settings);
         m_infoasscore.saveSettingsTo(settings);
         m_mafthreshold.saveSettingsTo(settings);
         m_memfree.saveSettingsTo(settings);
         m_memtotal.saveSettingsTo(settings);
         m_method.saveSettingsTo(settings);
         m_minqueryfrac.saveSettingsTo(settings);
         m_neargene.saveSettingsTo(settings);
         m_normscorethreshold.saveSettingsTo(settings);
         m_otherinfo.saveSettingsTo(settings);
         m_outfile.saveSettingsTo(settings);
         m_path2annovar.saveSettingsTo(settings);
         m_precedence.saveSettingsTo(settings);
         m_queryfile.saveSettingsTo(settings);
         m_rawscore.saveSettingsTo(settings);
         m_reverse.saveSettingsTo(settings);
         m_scorecolumn.saveSettingsTo(settings);
         m_scorethreshold.saveSettingsTo(settings);
         m_separate.saveSettingsTo(settings);
         m_seqpadding.saveSettingsTo(settings);
         m_siftthreshold.saveSettingsTo(settings);
         m_splicingthreshold.saveSettingsTo(settings);
         m_tablename.saveSettingsTo(settings);
         m_transcriptfunction.saveSettingsTo(settings);
         m_usechromosome.saveSettingsTo(settings);
         m_useindelsplicingthreshold.saveSettingsTo(settings);
         m_useneargene.saveSettingsTo(settings);
         m_usenormscorethreshold.saveSettingsTo(settings);
         m_usescorethreshold.saveSettingsTo(settings);
         m_usesplicingthreshold.saveSettingsTo(settings);
//         m_useoutfile.saveSettingsTo(settings);
//        m_usedbtype.saveSettingsTo(settings);
//         m_usebuildver.saveSettingsTo(settings);
         m_usegff3dbfile.saveSettingsTo(settings);
         m_usegenericdbfile.saveSettingsTo(settings);
         m_usevcfdbfile.saveSettingsTo(settings);
         m_usebedfile.saveSettingsTo(settings);
         m_usecolswanted.saveSettingsTo(settings);
         m_usescorecolumn.saveSettingsTo(settings);
         m_usetablename.saveSettingsTo(settings);
  //       m_usequeryfile.saveSettingsTo(settings);
         m_usesiftthreshold.saveSettingsTo(settings);
         m_vcfdbfile.saveSettingsTo(settings);
  //       m_webfrom.saveSettingsTo(settings);

         
         
     
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_batchsize.loadSettingsFrom(settings);
        m_bedfile.loadSettingsFrom(settings);
        m_buildver.loadSettingsFrom(settings);
        m_chromosome.loadSettingsFrom(settings);
        m_colswanted.loadSettingsFrom(settings);
        m_comment.loadSettingsFrom(settings);
        m_databaselocation.loadSettingsFrom(settings);
        m_dbtype.loadSettingsFrom(settings);
        m_exonsort.loadSettingsFrom(settings);
        m_expandbin.loadSettingsFrom(settings);
        m_genericdbfile.loadSettingsFrom(settings);
        m_genomebinsize.loadSettingsFrom(settings);
        m_gff3dbfile.loadSettingsFrom(settings);
        m_hgvs.loadSettingsFrom(settings);
        m_indelsplicingthreshold.loadSettingsFrom(settings);
        m_indexfilterthreshold.loadSettingsFrom(settings);
        m_infoasscore.loadSettingsFrom(settings);
        m_mafthreshold.loadSettingsFrom(settings);
        m_memfree.loadSettingsFrom(settings);
        m_memtotal.loadSettingsFrom(settings);
        m_method.loadSettingsFrom(settings);
        m_minqueryfrac.loadSettingsFrom(settings);
        m_neargene.loadSettingsFrom(settings);
        m_normscorethreshold.loadSettingsFrom(settings);
        m_otherinfo.loadSettingsFrom(settings);
        m_outfile.loadSettingsFrom(settings);
        m_path2annovar.loadSettingsFrom(settings);
        m_precedence.loadSettingsFrom(settings);
        m_queryfile.loadSettingsFrom(settings);
        m_rawscore.loadSettingsFrom(settings);
        m_reverse.loadSettingsFrom(settings);
        m_scorecolumn.loadSettingsFrom(settings);
        m_scorethreshold.loadSettingsFrom(settings);
        m_separate.loadSettingsFrom(settings);
        m_seqpadding.loadSettingsFrom(settings);
        m_siftthreshold.loadSettingsFrom(settings);
        m_splicingthreshold.loadSettingsFrom(settings);
        m_tablename.loadSettingsFrom(settings);
        m_transcriptfunction.loadSettingsFrom(settings);
        m_usechromosome.loadSettingsFrom(settings);
        m_useindelsplicingthreshold.loadSettingsFrom(settings);
        m_useneargene.loadSettingsFrom(settings);
        m_usenormscorethreshold.loadSettingsFrom(settings);
        m_usescorethreshold.loadSettingsFrom(settings);
        m_usesplicingthreshold.loadSettingsFrom(settings);
 //       m_useoutfile.loadSettingsFrom(settings);
 //       m_usedbtype.loadSettingsFrom(settings);
 //       m_usebuildver.loadSettingsFrom(settings);
        m_usegff3dbfile.loadSettingsFrom(settings);
        m_usegenericdbfile.loadSettingsFrom(settings);
        m_usevcfdbfile.loadSettingsFrom(settings);
        m_usebedfile.loadSettingsFrom(settings);
        m_usecolswanted.loadSettingsFrom(settings);
        m_usescorecolumn.loadSettingsFrom(settings);
        m_usetablename.loadSettingsFrom(settings);
 //       m_usequeryfile.loadSettingsFrom(settings);
        m_usesiftthreshold.loadSettingsFrom(settings);
        m_vcfdbfile.loadSettingsFrom(settings);
  //      m_webfrom.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_batchsize.validateSettings(settings);
        m_bedfile.validateSettings(settings);
        m_buildver.validateSettings(settings);
        m_chromosome.validateSettings(settings);
        m_colswanted.validateSettings(settings);
        m_comment.validateSettings(settings);
        m_databaselocation.validateSettings(settings);
        m_dbtype.validateSettings(settings);
        m_exonsort.validateSettings(settings);
        m_expandbin.validateSettings(settings);
        m_genericdbfile.validateSettings(settings);
        m_genomebinsize.validateSettings(settings);
        m_gff3dbfile.validateSettings(settings);
        m_hgvs.validateSettings(settings);
        m_indelsplicingthreshold.validateSettings(settings);
        m_indexfilterthreshold.validateSettings(settings);
        m_infoasscore.validateSettings(settings);
        m_mafthreshold.validateSettings(settings);
        m_memfree.validateSettings(settings);
        m_memtotal.validateSettings(settings);
        m_method.validateSettings(settings);
        m_minqueryfrac.validateSettings(settings);
        m_neargene.validateSettings(settings);
        m_normscorethreshold.validateSettings(settings);
        m_otherinfo.validateSettings(settings);
        m_outfile.validateSettings(settings);
        m_path2annovar.validateSettings(settings);
        m_precedence.validateSettings(settings);
        m_queryfile.validateSettings(settings);
        m_rawscore.validateSettings(settings);
        m_reverse.validateSettings(settings);
        m_scorecolumn.validateSettings(settings);
        m_scorethreshold.validateSettings(settings);
        m_separate.validateSettings(settings);
        m_seqpadding.validateSettings(settings);
        m_siftthreshold.validateSettings(settings);
        m_splicingthreshold.validateSettings(settings);
        m_tablename.validateSettings(settings);
        m_transcriptfunction.validateSettings(settings);
        m_usechromosome.validateSettings(settings);
        m_useindelsplicingthreshold.validateSettings(settings);
        m_useneargene.validateSettings(settings);
        m_usenormscorethreshold.validateSettings(settings);
        m_usescorethreshold.validateSettings(settings);
        m_usesplicingthreshold.validateSettings(settings);
//       m_useoutfile.validateSettings(settings);
//        m_usedbtype.validateSettings(settings);
 //       m_usebuildver.validateSettings(settings);
        m_usegff3dbfile.validateSettings(settings);
        m_usegenericdbfile.validateSettings(settings);
        m_usevcfdbfile.validateSettings(settings);
        m_usebedfile.validateSettings(settings);
        m_usecolswanted.validateSettings(settings);
        m_usescorecolumn.validateSettings(settings);
        m_usetablename.validateSettings(settings);
  //      m_usequeryfile.validateSettings(settings);
        m_usesiftthreshold.validateSettings(settings);
        m_vcfdbfile.validateSettings(settings);
 //       m_webfrom.validateSettings(settings);
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


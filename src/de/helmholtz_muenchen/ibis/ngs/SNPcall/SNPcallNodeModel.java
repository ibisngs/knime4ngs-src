package de.helmholtz_muenchen.ibis.ngs.SNPcall;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.ngs.QSub;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;


/**
 * This is the model implementation of SNPcall.
 * 
 */
public class SNPcallNodeModel extends NodeModel {
    
	/*Mpileup
	 * Input Options:

       -6           assume the quality is in the Illumina-1.3+ encoding
       -A           count anomalous read pairs
       -B           disable BAQ computation
       -b FILE      list of input BAM files [null]
       -C INT       parameter for adjusting mapQ; 0 to disable [0]
       -d INT       max per-BAM depth to avoid excessive memory usage [250]
       -E           extended BAQ for higher sensitivity but lower specificity
       -f FILE      faidx indexed reference sequence file [null]
       -G FILE      exclude read groups listed in FILE [null]
       -l FILE      list of positions (chr pos) or regions (BED) [null]
       -M INT       cap mapping quality at INT [60]
       -r STR       region in which pileup is generated [null]
       -R           ignore RG tags
       -q INT       skip alignments with mapQ smaller than INT [0]
       -Q INT       skip bases with baseQ/BAQ smaller than INT [13]

##SNP Call specific parameters
       -e INT       Phred-scaled gap extension seq error probability [20]
       -F FLOAT     minimum fraction of gapped reads for candidates [0.002]
       -h INT       coefficient for homopolymer errors [100]
       -I           do not perform indel calling
       -L INT       max per-sample depth for INDEL calling [250]
       -m INT       minimum gapped reads for indel candidates [1]
       -o INT       Phred-scaled gap open sequencing error probability [40]
       -P STR       comma separated list of platforms for indels [all]
            
       
 ##varFilter
Options: -Q INT    minimum RMS mapping quality for SNPs [10]
         -d INT    minimum read depth [2]
         -D INT    maximum read depth [10000000]
         -a INT    minimum number of alternate bases [2]
         -w INT    SNP within INT bp around a gap to be filtered [3]
         -W INT    window size for filtering adjacent gaps [10]
         -1 FLOAT  min P-value for strand bias (given PV4) [0.0001]
         -2 FLOAT  min P-value for baseQ bias [1e-100]
         -3 FLOAT  min P-value for mapQ bias [0]
         -4 FLOAT  min P-value for end distance bias [0.0001]
         -e FLOAT  min P-value for HWE (plus F<0) [0.0001]
         -p        print filtered variants
	 */

	//Mpileup
	public static final String CFGKEY_ENCODING = "encoding";
	public static final String CFGKEY_ANAMALOUS = "anamalous";
	public static final String CFGKEY_PROBREALIGN="probrealign";
	public static final String CFGKEY_DOWNGRADE ="downgrade";
	public static final String CFGKEY_MAXREADS="maxreads";
	public static final String CFGKEY_EXTENDBAQ="extendbaq";
//	public static final String CFGKEY_FAIDX="faidx";
	public static final String CFGKEY_BEDFILE="bedfile";
	public static final String CFGKEY_MINMAPQUAL="minmapqual";
	public static final String CFGKEY_MINBASEQUAL="minbasequal";
	public static final String CFGKEY_GAPEXTEND="gapextend";
	public static final String CFGKEY_HOMOPOLY="homopoly";
	public static final String CFGKEY_NOINDEL="noindel";
	public static final String CFGKEY_SKIPINDEL="skipindel";
	public static final String CFGKEY_GAPOPEN="gapopen";
	public static final String CFGKEY_IFBEDFILE="ifbedfile";	
	public static final String CFGKEY_EXCLUDEREADS="excludereads";
	public static final String CFGKEY_CAPMAPQUAL="capmapqual";
	public static final String CFGKEY_IGNORERG="ignorerg";
	public static final String CFGKEY_IFREADFILE="ifreadfile";	
	public static final String CFGKEY_MINFRAC="minfrac";
	public static final String CFGKEY_MINGAPREADS="mingapreads";
	
	public static final String CFGKEY_PILEUPREGION="pileuporegion";
	public static final String CFGKEY_PLATFORMLIST="platformlist";
	//varFilter
	public static final String CFGKEY_MINRMS="minrms";
	public static final String CFGKEY_MINREADDEPTH="minreaddepth";
	public static final String CFGKEY_MAXREADDEPTH="maxreaddepth";
	public static final String CFGKEY_MINALTBASE="minaltbases";
	public static final String CFGKEY_GAPFILTER="gapfilter";
	public static final String CFGKEY_ADJACENTGAPS="adjacentgaps";
	public static final String CFGKEY_STRANDPVAL="strandpval";
	public static final String CFGKEY_BASEQPVAL="basqpval";
	public static final String CFGKEY_MAPQPVAL="mapqpval";
	public static final String CFGKEY_ENDDISTPVAL="enddistpval";
	public static final String CFGKEY_HWEPVAL="hwepval";
	public static final String CFGKEY_PRINTFILTERED="printfiltered";
	
	//Mpileup
	public static final int DEFAULT_DOWNGRADE = 0;
	public static final int DEFAULT_MAXREADS = 250;
	public static final int DEFAULT_MINMAPQUAL=0;
	public static final int DEFAULT_MINBASEQUAL=13;
	public static final int DEFAULT_GAPEXTEND=20;
	public static final int DEFAULT_HOMOPOLY=100;
	public static final int DEFAULT_SKIPINDEL=250;
	public static final int DEFAULT_GAPOPEN=40;
	public static final int DEFAULT_CAPMAPQUAL=60;
	public static final double DEFAULT_MINFRAC=0.002;
	public static final int DEFAULT_MINGAPREADS=1;
    //varFilter
	public static final int DEFAULT_MINRMS=10;
	public static final int DEFAULT_MINREADDEPTH=2;
	public static final int DEFAULT_MAXREADDEPTH=10000000;
	public static final int DEFAULT_MINALTBASE=2;
	public static final int DEFAULT_GAPFILTER=3;
	public static final int DEFAULT_ADJACENTGAPS=10;
	public static final double DEFAULT_STRANDPVAL= 0.0001;
	public static final int DEFAULT_BASEQPVAL=100;
	public static final double DEFAULT_MAPQPVAL= 0;
	public static final double DEFAULT_ENDDISTPVAL= 0.0001;
	public static final double DEFAULT_HWEPVAL= 0.0001;
	
	
	
	
	//Mpileup
	private final SettingsModelBoolean m_encoding = new SettingsModelBoolean(
			CFGKEY_ENCODING, false);
	private final SettingsModelBoolean m_anamalous = new SettingsModelBoolean(
			CFGKEY_ANAMALOUS, false);
	private final SettingsModelBoolean m_probrealign = new SettingsModelBoolean(
			CFGKEY_PROBREALIGN, false);
	private final SettingsModelIntegerBounded m_downgrade = new SettingsModelIntegerBounded(
			CFGKEY_DOWNGRADE,DEFAULT_DOWNGRADE,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_maxreads = new SettingsModelIntegerBounded(
			CFGKEY_MAXREADS,DEFAULT_MAXREADS,0,Integer.MAX_VALUE);
	private final SettingsModelBoolean m_extendbaq = new SettingsModelBoolean(
			CFGKEY_EXTENDBAQ, false);
//	private final SettingsModelString m_faidx = new SettingsModelString(
//			SNPcallNodeModel.CFGEKY_FAIDX,"");
	private final SettingsModelString m_bedfile = new SettingsModelString(
			SNPcallNodeModel.CFGKEY_BEDFILE,"");
	private final SettingsModelIntegerBounded m_minmapqual = new SettingsModelIntegerBounded(
			CFGKEY_MINMAPQUAL,DEFAULT_MINMAPQUAL,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_minbasequal = new SettingsModelIntegerBounded(
			CFGKEY_MINBASEQUAL,DEFAULT_MINBASEQUAL,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_gapextend = new SettingsModelIntegerBounded(
			CFGKEY_GAPEXTEND,DEFAULT_GAPEXTEND,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_homopoly = new SettingsModelIntegerBounded(
			CFGKEY_HOMOPOLY,DEFAULT_HOMOPOLY,0,Integer.MAX_VALUE);
	private final SettingsModelBoolean m_noindel = new SettingsModelBoolean(
			CFGKEY_NOINDEL, false);
	private final SettingsModelIntegerBounded m_skipindel = new SettingsModelIntegerBounded(
			CFGKEY_SKIPINDEL,DEFAULT_SKIPINDEL,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_gapopen = new SettingsModelIntegerBounded(
			CFGKEY_GAPOPEN,DEFAULT_GAPOPEN,0,Integer.MAX_VALUE);
	//private final SettingsModelBoolean m_ifbedfile = new SettingsModelBoolean(
	//		CFGKEY_IFBEDFILE, false);
	private final SettingsModelString m_excludereadsfile = new SettingsModelString(
			SNPcallNodeModel.CFGKEY_EXCLUDEREADS,"");
	private final SettingsModelIntegerBounded m_capmapqual = new SettingsModelIntegerBounded(
			CFGKEY_CAPMAPQUAL,DEFAULT_CAPMAPQUAL,0,Integer.MAX_VALUE);
	private final SettingsModelBoolean m_ignorerg = new SettingsModelBoolean(
			CFGKEY_IGNORERG, false);
	private final SettingsModelDoubleBounded m_minfrac = new SettingsModelDoubleBounded(
			CFGKEY_MINFRAC,DEFAULT_MINFRAC,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_mingapreads = new SettingsModelIntegerBounded(
			CFGKEY_MINGAPREADS,DEFAULT_MINGAPREADS,0,Integer.MAX_VALUE);
	private final SettingsModelOptionalString m_pileupregion= new SettingsModelOptionalString(CFGKEY_PILEUPREGION, "", false);
	private final SettingsModelOptionalString m_platformlist= new SettingsModelOptionalString(CFGKEY_PLATFORMLIST, "", false);
	//varFilter	
	private final SettingsModelIntegerBounded m_minrms = new SettingsModelIntegerBounded(
			CFGKEY_MINRMS,DEFAULT_MINRMS,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_minreaddepth = new SettingsModelIntegerBounded(
			CFGKEY_MINREADDEPTH,DEFAULT_MINREADDEPTH,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_maxreaddepth = new SettingsModelIntegerBounded(
			CFGKEY_MAXREADDEPTH,DEFAULT_MAXREADDEPTH,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_minaltbase = new SettingsModelIntegerBounded(
			CFGKEY_MINALTBASE,DEFAULT_MINALTBASE,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_gapfilter = new SettingsModelIntegerBounded(
			CFGKEY_GAPFILTER,DEFAULT_GAPFILTER,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_adjacentgaps = new SettingsModelIntegerBounded(
			CFGKEY_ADJACENTGAPS,DEFAULT_ADJACENTGAPS,0,Integer.MAX_VALUE);
	private final SettingsModelDoubleBounded m_strandpval = new SettingsModelDoubleBounded(
			CFGKEY_STRANDPVAL,DEFAULT_STRANDPVAL,0,Double.MAX_VALUE);
	private final SettingsModelIntegerBounded m_baseqpval = new SettingsModelIntegerBounded(
			CFGKEY_BASEQPVAL,DEFAULT_BASEQPVAL,0,Integer.MAX_VALUE);
	private final SettingsModelDoubleBounded m_mapqpval = new SettingsModelDoubleBounded(
			CFGKEY_MAPQPVAL,DEFAULT_MAPQPVAL,0,Double.MAX_VALUE);
	private final SettingsModelDoubleBounded m_enddistpval = new SettingsModelDoubleBounded(
			CFGKEY_ENDDISTPVAL,DEFAULT_ENDDISTPVAL,0,Double.MAX_VALUE);
	private final SettingsModelDoubleBounded m_hwepval = new SettingsModelDoubleBounded(
			CFGKEY_HWEPVAL,DEFAULT_HWEPVAL,0,Double.MAX_VALUE);	
	private final SettingsModelBoolean m_printfiltered = new SettingsModelBoolean(
					CFGKEY_PRINTFILTERED, false);
	private final SettingsModelBoolean m_ifbedfile = new SettingsModelBoolean(
			CFGKEY_IFBEDFILE, false);
	private final SettingsModelBoolean m_ifreadfile = new SettingsModelBoolean(
			CFGKEY_IFREADFILE, false);
	

    /**
     * Constructor for the node model.
     */
    protected SNPcallNodeModel() {
    	
        super(1, 0);
        
    	m_bedfile.setEnabled(false);
    	m_excludereadsfile.setEnabled(false);
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String path2bamfile = inData[0].iterator().next().getCell(1).toString();

    	/**Initialize logfile**/
    	String logfile = path2bamfile.substring(0,path2bamfile.lastIndexOf("/")+1)+"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("SNPcall"));
    	/**end initializing logfile**/
    	
    	String path2samtools = inData[0].iterator().next().getCell(0).toString();
    	String path2bcftools = path2samtools.substring(0,path2samtools.length()-8)+"bcftools/bcftools";
    	String path2vcfutils = path2samtools.substring(0,path2samtools.length()-8)+"bcftools/vcfutils.pl";
    	String path2seqfile;
    	if(inData[0].iterator().next().getNumCells() > 2) {
    		path2seqfile = inData[0].iterator().next().getCell(2).toString();
    	} else {
    		path2seqfile = getAvailableInputFlowVariables().get("Path2seqFile").getStringValue();
    	}
    	String baseName = path2bamfile;
    	String path2outputfile = path2bamfile.substring(0,path2bamfile.length()-4)+"_SNPs_raw.bcf";
    	String path2outputfile2 = path2bamfile.substring(0,path2bamfile.length()-4)+"_SNPs_flt.vcf";
    	
    	String fadix = " faidx ";
    	String mileup = " mpileup ";
    	String view = " view -bvcg ";
    	//String varfilter = " varFilter -D100 > ";
    	String varfilter = " varFilter ";
    	
    /**
     * 	Process One:
     * Index reference sequence file
     */
    	ShowOutput.writeLogFile("Index reference sequence");
    	String com = path2samtools + fadix + path2seqfile;
    	// begin QueueSub #################################################
		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_SNPcall-faidx";
			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
			String logfle = path2bamfile.substring(0,path2bamfile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
			new QSub(com, name, memory, logfle, true);
 			logBuffer.append("QSub: " + com + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
		// end QueueSub ###################################################
		} else {
	    	System.out.println(com);
	    	Process p_fadix = Runtime.getRuntime().exec(com);
	    	p_fadix.waitFor();
	    	logBuffer.append(ShowOutput.getLogEntry(p_fadix, com));
		}
    	
  	/**
  	 * Process Two:
  	 * Run mpileup 	
  	 */
   // samtools mpileup -uf ref.fa aln1.bam aln2.bam | bcftools view -bvcg - > var.raw.bcf
    	System.out.println("Mpileup Process and Bcftools");
    	String encoding="";
    	if(m_encoding.getBooleanValue()){encoding="-6 ";}
    	String anamalous="";
    	if(m_anamalous.getBooleanValue()){anamalous="-A ";}
    	String probrealgin="";
    	if(m_probrealign.getBooleanValue()){probrealgin="-B ";}
    	String downgrade = "-C "+m_downgrade.getIntValue()+" ";
    	String maxreads = "-d "+m_maxreads.getIntValue()+" ";
    	String extendbaq="";
    	if(m_extendbaq.getBooleanValue()){extendbaq="-E ";}
    	String faidx2 = "-f "+path2seqfile+" ";
    	String bedfile="";
    	if(m_bedfile.isEnabled()){bedfile = "-l "+m_bedfile.getStringValue()+" ";}
    	String minmapqual = "-q "+m_minmapqual.getIntValue()+" ";
    	String minbasequal = "-Q "+m_minbasequal.getIntValue()+" ";
    	String output ="-g ";
    	String gapextend="";
    	if(m_gapextend.isEnabled()){gapextend = "-e "+m_gapextend.getIntValue()+" ";}
    	String minfrac="";
    	if(m_minfrac.isEnabled()){minfrac="-F "+m_minfrac.getDoubleValue()+" ";}
    	String homopoly="";
    	if(m_homopoly.isEnabled()){homopoly = "-h "+m_homopoly.getIntValue()+" ";}
    	String noindel="";
    	if(m_noindel.isEnabled() && m_noindel.getBooleanValue()){noindel="-I ";}
    	String skipindel ="";
    	if(m_skipindel.isEnabled()){skipindel = "-L "+m_skipindel.getIntValue()+" ";}
    	String mingapreads="";
    	if(m_mingapreads.isEnabled()){mingapreads="-m "+m_mingapreads.getIntValue()+" ";}
    	String gapopen="";
    	if(m_gapopen.isEnabled()){gapopen = "-o "+m_gapopen.getIntValue()+" ";}
    	String excludereads="";
    	if(m_excludereadsfile.isEnabled()){excludereads="-G "+m_excludereadsfile.getStringValue()+" ";}
    	String capmapqual="-M "+m_capmapqual.getIntValue()+" ";
    	String ignorerg="";
    	if(m_ignorerg.getBooleanValue()){ignorerg="-R ";}
    	String pileupregion="";
    	if(m_pileupregion.isActive()){pileupregion="-r "+m_pileupregion.getStringValue()+" ";}
    	String platformlist="";
    	if(m_platformlist.isActive()){platformlist="-P "+m_platformlist.getStringValue()+" ";}
    	String mptempfile = baseName + "_mpileup.bcf";
	
    	com = path2samtools + mileup +encoding+anamalous+probrealgin+downgrade+maxreads+extendbaq+excludereads+bedfile+minmapqual+capmapqual+ignorerg+pileupregion+minbasequal+output+gapextend+minfrac+homopoly+noindel+skipindel+mingapreads+gapopen+platformlist+faidx2+path2bamfile + " > " + mptempfile;
    	// begin QueueSub #################################################
		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_SNPcall-mpileup";
			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
			String logfle = path2bamfile.substring(0,path2bamfile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
			new QSub(com, name, memory, logfle, true);
 			logBuffer.append("QSub: " + com + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
		// end QueueSub ###################################################
		} else {
	    	System.out.println(com);
	    	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", com);
	    	Process p_mileup = b.start();
	    	p_mileup.waitFor();
	    	logBuffer.append(ShowOutput.getLogEntry(p_mileup, com));
		}
    	
    	com = path2bcftools + view + mptempfile + " > " + path2outputfile;
    	// begin QueueSub #################################################
		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_SNPcall-view";
			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
			String logfle = path2bamfile.substring(0,path2bamfile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
			new QSub(com, name, memory, logfle, true);
 			logBuffer.append("QSub: " + com + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
		// end QueueSub ###################################################
		} else {
	    	System.out.println(com);
	    	ProcessBuilder b1 = new ProcessBuilder("/bin/sh", "-c", com);
	    	Process p_bcfTview = b1.start();
	    	p_bcfTview.waitFor();
	    	logBuffer.append(ShowOutput.getLogEntry(p_bcfTview, com));
		}
    	
    /**
     * Process Three:
     * Run bcftools and varFilter	
     */    	
    // bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf
    	String minrms="-Q "+m_minrms.getIntValue()+" ";
    	String minreaddepth="-d "+m_minreaddepth.getIntValue()+" ";
    	String maxreaddepth="-D "+m_maxreaddepth.getIntValue()+" ";
    	String minaltbases="-a "+m_minaltbase.getIntValue()+" ";
    	String gapfilter="-w "+m_gapfilter.getIntValue()+" ";
    	String adjacentgaps="-W "+m_adjacentgaps.getIntValue()+" ";
    	String strandpval="-1 "+m_strandpval.getDoubleValue()+" ";
    	String baseqpval = "-2 1.0E-" + m_baseqpval.getIntValue() + " ";
    	String mapqpval="-3 "+m_mapqpval.getDoubleValue()+" ";
    	String enddistpval = "-4 "+m_enddistpval.getDoubleValue()+" ";
    	String hwepval="-e "+m_hwepval.getDoubleValue()+" ";
    	String printfiltered="";
    	if(m_printfiltered.getBooleanValue()){printfiltered="-p ";}	
    	String outvcf = path2outputfile + ".vcf";
    
    	
    	logBuffer.append("Bcftools view and Vcfutils");
    	com =  path2bcftools +" view "+ path2outputfile + " > " + outvcf;
    	// begin QueueSub #################################################
		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_SNPcall-view";
			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
			String logfle = path2bamfile.substring(0,path2bamfile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
			new QSub(com, name, memory, logfle, true);
 			logBuffer.append("QSub: " + com + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
		// end QueueSub ###################################################
		} else {
	    	System.out.println(com);
	    	ProcessBuilder b2 = new ProcessBuilder("/bin/sh", "-c", com);
	    	Process p_filter = b2.start();
	    	p_filter.waitFor();
	    	logBuffer.append(ShowOutput.getLogEntry(p_filter, com));
		}
    	
    	com = path2vcfutils + varfilter +minrms+minreaddepth+maxreaddepth+minaltbases+gapfilter+adjacentgaps+strandpval+baseqpval+mapqpval+enddistpval+hwepval+printfiltered+outvcf+" > " + path2outputfile2;
    	// begin QueueSub #################################################
		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_SNPcall-varFilter";
			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
			String logfle = path2bamfile.substring(0,path2bamfile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
			new QSub(com, name, memory, logfle, true);
 			logBuffer.append("QSub: " + com + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
		// end QueueSub ###################################################
		} else {
	    	System.out.println(com);
	    	ProcessBuilder b2 = new ProcessBuilder("/bin/sh", "-c", com);
	    	Process p_filter = b2.start();
	    	p_filter.waitFor();
	    	logBuffer.append(ShowOutput.getLogEntry(p_filter, com));
		}
    	
    	logBuffer.append(ShowOutput.getNodeEndTime());
    	ShowOutput.writeLogFile(logBuffer);
    	
        return new BufferedDataTable[]{};
    }

    /*private void deleteTemporaryFile(String filename) {
    	File file = new File(filename);    	        
		if(file.exists()){
			file.delete();
		}
    }*/
    
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
    	
    	// Check input ports
    	String[] cn=inSpecs[0].getColumnNames();
    	if(!cn[0].equals("") && !cn[1].equals("") && !cn[0].equals("Path2SamTools") && !cn[1].equals("Path2BAMFile")) {
    		throw new InvalidSettingsException("This node is incompatible with the previous node. The outport of the previous node has to fit to the inport of this node.");
    	}
    	
        return new DataTableSpec[]{};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_anamalous.saveSettingsTo(settings);
    	m_bedfile.saveSettingsTo(settings);
    	m_downgrade.saveSettingsTo(settings);
    	m_encoding.saveSettingsTo(settings);
    	m_extendbaq.saveSettingsTo(settings);
    	m_maxreads.saveSettingsTo(settings);
    	m_minbasequal.saveSettingsTo(settings);
    	m_minmapqual.saveSettingsTo(settings);
    	m_probrealign.saveSettingsTo(settings);
    	m_gapextend.saveSettingsTo(settings);
    	m_homopoly.saveSettingsTo(settings);
    	m_noindel.saveSettingsTo(settings);
    	m_skipindel.saveSettingsTo(settings);
    	m_gapopen.saveSettingsTo(settings);
    	m_ifbedfile.saveSettingsTo(settings);
    	m_ifreadfile.saveSettingsTo(settings);
    	m_excludereadsfile.saveSettingsTo(settings);
    	m_capmapqual.saveSettingsTo(settings);
    	m_ignorerg.saveSettingsTo(settings);
    	m_minfrac.saveSettingsTo(settings);
    	m_mingapreads.saveSettingsTo(settings);
    	m_pileupregion.saveSettingsTo(settings);
    	m_platformlist.saveSettingsTo(settings);
    	
    	m_minrms.saveSettingsTo(settings);
    	m_minreaddepth.saveSettingsTo(settings);
    	m_maxreaddepth.saveSettingsTo(settings);
    	m_minaltbase.saveSettingsTo(settings);
    	m_gapfilter.saveSettingsTo(settings);
    	m_adjacentgaps.saveSettingsTo(settings);
    	m_strandpval.saveSettingsTo(settings);
    	m_baseqpval.saveSettingsTo(settings);
    	m_mapqpval.saveSettingsTo(settings);
    	m_enddistpval.saveSettingsTo(settings);
    	m_hwepval.saveSettingsTo(settings);
    	m_printfiltered.saveSettingsTo(settings);
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        	m_anamalous.loadSettingsFrom(settings);
        	m_bedfile.loadSettingsFrom(settings);
        	m_downgrade.loadSettingsFrom(settings);
        	m_encoding.loadSettingsFrom(settings);
        	m_extendbaq.loadSettingsFrom(settings);
        	m_maxreads.loadSettingsFrom(settings);
        	m_minbasequal.loadSettingsFrom(settings);
        	m_minmapqual.loadSettingsFrom(settings);
        	m_probrealign.loadSettingsFrom(settings);
        	m_gapextend.loadSettingsFrom(settings);
        	m_homopoly.loadSettingsFrom(settings);
        	m_noindel.loadSettingsFrom(settings);
        	m_skipindel.loadSettingsFrom(settings);
        	m_gapopen.loadSettingsFrom(settings);
        	m_ifbedfile.loadSettingsFrom(settings);
        	m_ifreadfile.loadSettingsFrom(settings);
        	m_excludereadsfile.loadSettingsFrom(settings);
        	m_capmapqual.loadSettingsFrom(settings);
        	m_ignorerg.loadSettingsFrom(settings);
        	m_minfrac.loadSettingsFrom(settings);
        	m_mingapreads.loadSettingsFrom(settings);
        	m_pileupregion.loadSettingsFrom(settings);
        	m_platformlist.loadSettingsFrom(settings);
        	
        	m_minrms.loadSettingsFrom(settings);
        	m_minreaddepth.loadSettingsFrom(settings);
        	m_maxreaddepth.loadSettingsFrom(settings);
        	m_minaltbase.loadSettingsFrom(settings);
        	m_gapfilter.loadSettingsFrom(settings);
        	m_adjacentgaps.loadSettingsFrom(settings);
        	m_strandpval.loadSettingsFrom(settings);
        	m_baseqpval.loadSettingsFrom(settings);
        	m_mapqpval.loadSettingsFrom(settings);
        	m_enddistpval.loadSettingsFrom(settings);
        	m_hwepval.loadSettingsFrom(settings);
        	m_printfiltered.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        	m_anamalous.validateSettings(settings);
        	m_bedfile.validateSettings(settings);
        	m_downgrade.validateSettings(settings);
        	m_encoding.validateSettings(settings);
        	m_extendbaq.validateSettings(settings);
        	m_maxreads.validateSettings(settings);
        	m_minbasequal.validateSettings(settings);
        	m_minmapqual.validateSettings(settings);
        	m_probrealign.validateSettings(settings);
        	m_gapextend.validateSettings(settings);
        	m_homopoly.validateSettings(settings);
        	m_noindel.validateSettings(settings);
        	m_skipindel.validateSettings(settings);
        	m_gapopen.validateSettings(settings);
        	m_ifbedfile.validateSettings(settings);
        	m_ifreadfile.validateSettings(settings);
        	m_excludereadsfile.validateSettings(settings);
        	m_capmapqual.validateSettings(settings);
        	m_ignorerg.validateSettings(settings);
        	m_minfrac.validateSettings(settings);
        	m_mingapreads.validateSettings(settings);
        	m_pileupregion.validateSettings(settings);
        	m_platformlist.validateSettings(settings);
        	
        	m_minrms.validateSettings(settings);
        	m_minreaddepth.validateSettings(settings);
        	m_maxreaddepth.validateSettings(settings);
        	m_minaltbase.validateSettings(settings);
        	m_gapfilter.validateSettings(settings);
        	m_adjacentgaps.validateSettings(settings);
        	m_strandpval.validateSettings(settings);
        	m_baseqpval.validateSettings(settings);
        	m_mapqpval.validateSettings(settings);
        	m_enddistpval.validateSettings(settings);
        	m_hwepval.validateSettings(settings);
        	m_printfiltered.validateSettings(settings);
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


package de.helmholtz_muenchen.ibis.ngs.mpileup;

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
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of Mpileup.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class MpileupNodeModel extends NodeModel {
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

Output options:

       -D           output per-sample DP in BCF (require -g/-u)
       -g           generate BCF output (genotype likelihoods)
       -O           output base positions on reads (disabled by -g/-u)
       -s           output mapping quality (disabled by -g/-u)
       -S           output per-sample strand bias P-value in BCF (require -g/-u)
       -u           generate uncompress BCF output



##SNP Call specific parameters
       -e INT       Phred-scaled gap extension seq error probability [20]
       -F FLOAT     minimum fraction of gapped reads for candidates [0.002]
       -h INT       coefficient for homopolymer errors [100]
       -I           do not perform indel calling
       -L INT       max per-sample depth for INDEL calling [250]
       -m INT       minimum gapped reads for indel candidates [1]
       -o INT       Phred-scaled gap open sequencing error probability [40]
       -P STR       comma separated list of platforms for indels [all]			

	 */

	//Mpileup
	public static final String CFGKEY_ENCODING = "encoding";
	public static final String CFGKEY_ANAMALOUS = "anamalous";
	public static final String CFGKEY_PROBREALIGN="probrealign";
	public static final String CFGKEY_DOWNGRADE ="downgrade";
	public static final String CFGKEY_MAXREADS="maxreads";
	public static final String CFGKEY_EXTENDBAQ="extendbaq";
	public static final String CFGKEY_USEFAIDX="usefaidx";
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
	
	public static final String CFGKEY_PILEUPREGION="pileuoregion";
	public static final String CFGKEY_PLATFORMLIST="platformlist";
	
	public static final String CFGKEY_OUTPERSAMPLE="outpersample";
	public static final String CFGKEY_BCFOUTPUT="bcfoutput";
	public static final String CFGKEY_OUTBASEPOSITIONS="basepositions";
	public static final String CFGKEY_OUTMAPQUAL="outmapqual";
	public static final String CFGKEY_STRANDBIASPVAL="strandbiaspval";
	public static final String CFGKEY_UNCOMPRESSEDBCF="uncompressedbcf";
	
	public static final String CFGKEY_IFINDEX="ifindex";
	
	
	//Mpileup
	public static final int DEFAULT_DOWNGRADE = 0;
	public static final int DEFAULT_MAXREADS = 8000;
	public static final int DEFAULT_MINMAPQUAL=0;
	public static final int DEFAULT_MINBASEQUAL=13;
	public static final int DEFAULT_GAPEXTEND=20;
	public static final int DEFAULT_HOMOPOLY=100;
	public static final int DEFAULT_SKIPINDEL=250;
	public static final int DEFAULT_GAPOPEN=40;
	public static final int DEFAULT_CAPMAPQUAL=60;
	public static final double DEFAULT_MINFRAC=0.002;
	public static final int DEFAULT_MINGAPREADS=1;

	
	/**
	 * Input models
	 */
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
	private final SettingsModelString m_excludereadsfile = new SettingsModelString(
			MpileupNodeModel.CFGKEY_EXCLUDEREADS,"");
	private final SettingsModelString m_bedfile = new SettingsModelString(
			MpileupNodeModel.CFGKEY_BEDFILE,"");
	private final SettingsModelIntegerBounded m_capmapqual = new SettingsModelIntegerBounded(
			CFGKEY_CAPMAPQUAL,DEFAULT_CAPMAPQUAL,0,Integer.MAX_VALUE);
	private final SettingsModelBoolean m_ignorerg = new SettingsModelBoolean(
			CFGKEY_IGNORERG, false);
	private final SettingsModelIntegerBounded m_minmapqual = new SettingsModelIntegerBounded(
			CFGKEY_MINMAPQUAL,DEFAULT_MINMAPQUAL,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_minbasequal = new SettingsModelIntegerBounded(
			CFGKEY_MINBASEQUAL,DEFAULT_MINBASEQUAL,0,Integer.MAX_VALUE);
	private final SettingsModelBoolean m_ifreadsfile = new SettingsModelBoolean(
			MpileupNodeModel.CFGKEY_IFREADFILE,false);
	private final SettingsModelBoolean m_ifbedfile = new SettingsModelBoolean(
			MpileupNodeModel.CFGKEY_IFBEDFILE,false);
	private final SettingsModelBoolean m_usefaidxfile = new SettingsModelBoolean(
			MpileupNodeModel.CFGKEY_USEFAIDX,true);
	
	private final SettingsModelOptionalString m_pileupregion= new SettingsModelOptionalString(CFGKEY_PILEUPREGION, "", false);

	
	
	/**
	 *SNP Call specific parameters
	 */
	private final SettingsModelIntegerBounded m_gapextend = new SettingsModelIntegerBounded(
			CFGKEY_GAPEXTEND,DEFAULT_GAPEXTEND,0,Integer.MAX_VALUE);
	private final SettingsModelDoubleBounded m_minfrac = new SettingsModelDoubleBounded(
			CFGKEY_MINFRAC,DEFAULT_MINFRAC,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_homopoly = new SettingsModelIntegerBounded(
			CFGKEY_HOMOPOLY,DEFAULT_HOMOPOLY,0,Integer.MAX_VALUE);
	private final SettingsModelBoolean m_noindel = new SettingsModelBoolean(
			CFGKEY_NOINDEL, false);
	private final SettingsModelIntegerBounded m_skipindel = new SettingsModelIntegerBounded(
			CFGKEY_SKIPINDEL,DEFAULT_SKIPINDEL,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_gapopen = new SettingsModelIntegerBounded(
			CFGKEY_GAPOPEN,DEFAULT_GAPOPEN,0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_mingapreads = new SettingsModelIntegerBounded(
			CFGKEY_MINGAPREADS,DEFAULT_MINGAPREADS,0,Integer.MAX_VALUE);
	private final SettingsModelOptionalString m_platformlist= new SettingsModelOptionalString(CFGKEY_PLATFORMLIST, "", false);
	
	
	/**
	 * Output models
	 */
	private final SettingsModelBoolean m_outpersample = new SettingsModelBoolean(
			MpileupNodeModel.CFGKEY_OUTPERSAMPLE,false);
	private final SettingsModelBoolean m_bcfoutput = new SettingsModelBoolean(
			MpileupNodeModel.CFGKEY_BCFOUTPUT,true);
	private final SettingsModelBoolean m_outbasepositions = new SettingsModelBoolean(
			MpileupNodeModel.CFGKEY_OUTBASEPOSITIONS,false);
	private final SettingsModelBoolean m_outmapqual = new SettingsModelBoolean(
			MpileupNodeModel.CFGKEY_OUTMAPQUAL,false);
	private final SettingsModelBoolean m_strandbiaspval = new SettingsModelBoolean(
			MpileupNodeModel.CFGKEY_STRANDBIASPVAL,false);
	private final SettingsModelBoolean m_uncompressedbcf = new SettingsModelBoolean(
			MpileupNodeModel.CFGKEY_UNCOMPRESSEDBCF,false);
	
	/**
	 * Extra models
	 */
	private final SettingsModelBoolean m_ifindex = new SettingsModelBoolean(
			MpileupNodeModel.CFGKEY_IFINDEX,true);
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(MpileupNodeModel.class);
	
	//The Output Col Names
	public static final String OUT_COL1 = "Path2SamTools";
	public static final String OUT_COL2 = "Path2MpileupOutfile";
	
    /**
     * Constructor for the node model.
     */
    protected MpileupNodeModel() {
    	
        super(1, 1);
        
    	m_bedfile.setEnabled(false);
    	m_excludereadsfile.setEnabled(false);
    	m_outmapqual.setEnabled(false);
    	m_outbasepositions.setEnabled(false);
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    		String path2seqfile = getAvailableInputFlowVariables().get("Path2seqFile").getStringValue();
        	String path2samtools = inData[0].iterator().next().getCell(0).toString();
        	String path2bamfile = inData[0].iterator().next().getCell(1).toString();
    	
	    	
	    	/**Initialize logfile**/
	    	String logfile = path2bamfile.substring(0,path2bamfile.lastIndexOf("/")+1)+"logfile.txt";
	    	System.out.println(logfile);
	    	ShowOutput.setLogFile(logfile);
	    	StringBuffer logBuffer = new StringBuffer(50);
	    	logBuffer.append(ShowOutput.getNodeStartTime("Mpileup"));
	    	/**end initializing logfile**/

    	
	    	String outfile = "";
	    	
	    	//If indexing is needed
	    	if(m_ifindex.getBooleanValue()){
	    	    /**
	    	     * 	Process One:
	    	     * Index reference sequence file
	    	     */
	    			ArrayList<String> command = new ArrayList<String>();
	    	    	LOGGER.info("Index reference sequence");
	    	    	command.add(path2samtools + " faidx");
	    	    	command.add(path2seqfile);
	    	    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
	    			
	    	}

	    	
	  	/**
	  	 * Process Two:
	  	 * Run mpileup 	
	  	 */
		   // samtools mpileup -uf ref.fa aln1.bam aln2.bam
	    		ArrayList<String> command = new ArrayList<String>();
	    		LOGGER.info("Mpileup Process");
	    		command.add(path2samtools+" mpileup");
	    		command.add(path2bamfile);
	    		
		    	if(m_encoding.getBooleanValue()){command.add("-6");}
		    	if(m_anamalous.getBooleanValue()){command.add("-A");}
		    	if(m_probrealign.getBooleanValue()){command.add("-B");}
		    	command.add("-C "+m_downgrade.getIntValue());
		    	command.add("-d "+m_maxreads.getIntValue());
		    	if(m_extendbaq.getBooleanValue()){command.add("-E");}
		    	if(m_usefaidxfile.getBooleanValue()){command.add("-f "+path2seqfile);}
		    	if(m_bedfile.isEnabled()){command.add("-l "+m_bedfile.getStringValue());}
		    	command.add("-q "+m_minmapqual.getIntValue());
		    	command.add("-Q "+m_minbasequal.getIntValue());
		    	command.add("-M "+m_capmapqual.getIntValue());
		    	if(m_excludereadsfile.isEnabled()){command.add("-G "+m_excludereadsfile.getStringValue());}
		    	if(m_ignorerg.getBooleanValue()){command.add("-R");}
		    	if(m_pileupregion.isActive()){command.add("-r "+m_pileupregion.getStringValue());}
		
		    	
		    	if(m_outpersample.isEnabled() && m_outpersample.getBooleanValue()){command.add("-D");}
		    	if(m_bcfoutput.isEnabled() && m_bcfoutput.getBooleanValue()){command.add("-g");}
		    	if(m_outbasepositions.isEnabled() && m_outbasepositions.getBooleanValue()){command.add("-O");}
		    	if(m_outmapqual.isEnabled() && m_outmapqual.getBooleanValue()){command.add("-s");}
		    	if(m_strandbiaspval.isEnabled() && m_strandbiaspval.getBooleanValue()){command.add("-S");}
		    	if(m_uncompressedbcf.isEnabled() && m_uncompressedbcf.getBooleanValue()){command.add("-u");}
		    		
		    	if(m_gapextend.isEnabled()){command.add("-e "+m_gapextend.getIntValue());}
		    	if(m_minfrac.isEnabled()){command.add("-F "+m_minfrac.getDoubleValue());}
		    	if(m_homopoly.isEnabled()){command.add("-h "+m_homopoly.getIntValue());}
		    	if(m_noindel.isEnabled() && m_noindel.getBooleanValue()){command.add("-I");}
		    	if(m_skipindel.isEnabled()){command.add("-L "+m_skipindel.getIntValue());}
		    	if(m_mingapreads.isEnabled()){command.add("-m "+m_mingapreads.getIntValue());}
		    	if(m_gapopen.isEnabled()){command.add("-o "+m_gapopen.getIntValue());}
		    	if(m_platformlist.isActive()){command.add("-P "+m_platformlist.getStringValue());}
			
		    	//If bcf output
		    	if((m_bcfoutput.isEnabled() && m_bcfoutput.getBooleanValue()) || (m_uncompressedbcf.isEnabled() && m_uncompressedbcf.getBooleanValue())){
		    		outfile =path2bamfile+"_pileup.bcf";
		    	}else{
		       		outfile =path2bamfile+".pileup";
		    	}
		    	
		    	/**
		    	 * Execute
		    	 */
		    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER,outfile);
		    	logBuffer.append(ShowOutput.getNodeEndTime());
		    	ShowOutput.writeLogFile(logBuffer);
		    	
		    	/**
		    	 * Output
		    	 */
		    	BufferedDataContainer cont = exec.createDataContainer(
		    			new DataTableSpec(
		    			new DataColumnSpec[]{
		    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
		    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
		    	
		    	FileCell[] c = new FileCell[]{
		    			(FileCell) FileCellFactory.create(path2samtools),
		    			(FileCell) FileCellFactory.create(outfile)};
		    	
		    	cont.addRowToTable(new DefaultRow("Row0",c));
		    	cont.close();
		    	BufferedDataTable outTable = cont.getTable();
		    	

				return new BufferedDataTable[]{outTable};
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
    	
		/**
		 * Check if reference file is available
		 */
    	if(getAvailableInputFlowVariables().containsKey("Path2seqFile")){
    		String referencefile=getAvailableInputFlowVariables().get("Path2seqFile").getStringValue();
    		System.out.println(referencefile);
    		if(referencefile.equals("No reference file")){
    			setWarningMessage("Reference sequence file variable not available. Probably you used BAMSAMConverter as a stand-alone node. You cannot use reference sequence for Mpileup now.\n Please use BAMLoader if you want to specify a reference sequence.");
    			//No reference, set Dialog Options to disabled
    			System.out.println("No-->Reference--disable");
    			m_usefaidxfile.setEnabled(false);
    			m_usefaidxfile.setBooleanValue(false);
    			m_ifindex.setBooleanValue(false);
    			m_ifindex.setEnabled(false);
    		}else{
        		//Reference file available-->Enable dialog
    			m_usefaidxfile.setEnabled(true);
    			m_usefaidxfile.setBooleanValue(true);
    			m_ifindex.setEnabled(true);
    		}

    	}else{
    		throw new InvalidSettingsException("Reference file variable not available. Mpileup should be connected to BAMSAMConverter or BAMLoader.");
    	}
    	
    /**
     * Check In and Outfiles
     */
    	String outfile="";
    	String path2bamfile="";
    	if(getAvailableInputFlowVariables().containsKey("MPILEUPINFILE")){
    		//Get outfile from BAMLoader/BAMSAMConverter
        	path2bamfile=getAvailableInputFlowVariables().get("MPILEUPINFILE").getStringValue();
        	String suffix = path2bamfile.substring(path2bamfile.lastIndexOf("."));
        	if(!suffix.equals(".bam")){
        		throw new InvalidSettingsException("Something wrong. Infile not in .bam format. Instead it is in "+suffix+" format");
        	}
        	//Infile is bam...prepare output file
        	if((m_bcfoutput.isEnabled() && m_bcfoutput.getBooleanValue()) || (m_uncompressedbcf.isEnabled() && m_uncompressedbcf.getBooleanValue())){
        		
        		outfile =path2bamfile+"_pileup.bcf";
        	}else{
           		outfile =path2bamfile+".pileup";
        	}
    		/**
    		 * Check if Outfile exists
    		 */		
            	File outpath = new File(outfile);
            	if(outpath.exists()){
            		throw new InvalidSettingsException("Outfile "+outpath+" already exists ! Please rename or move to other directory.");
            	}
    	}else{
    		throw new InvalidSettingsException("Infile variable not available. Mpileup should be connected to BAMSAMConverter or BAMLoader.");
    	}
    	
        return new DataTableSpec[]{null};
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
    	m_ifreadsfile.saveSettingsTo(settings);
    	m_excludereadsfile.saveSettingsTo(settings);
    	m_capmapqual.saveSettingsTo(settings);
    	m_ignorerg.saveSettingsTo(settings);
    	m_minfrac.saveSettingsTo(settings);
    	m_mingapreads.saveSettingsTo(settings);
    	m_pileupregion.saveSettingsTo(settings);
    	m_platformlist.saveSettingsTo(settings);
    	
    	m_outpersample.saveSettingsTo(settings);
    	m_bcfoutput.saveSettingsTo(settings);
    	m_outbasepositions.saveSettingsTo(settings);
    	m_outmapqual.saveSettingsTo(settings);
    	m_strandbiaspval.saveSettingsTo(settings);
    	m_uncompressedbcf.saveSettingsTo(settings);
    	m_ifindex.saveSettingsTo(settings);
    	m_usefaidxfile.saveSettingsTo(settings);
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
    	m_ifreadsfile.loadSettingsFrom(settings);
    	m_excludereadsfile.loadSettingsFrom(settings);
    	m_capmapqual.loadSettingsFrom(settings);
    	m_ignorerg.loadSettingsFrom(settings);
    	m_minfrac.loadSettingsFrom(settings);
    	m_mingapreads.loadSettingsFrom(settings);
    	m_pileupregion.loadSettingsFrom(settings);
    	m_platformlist.loadSettingsFrom(settings);
    	
    	m_outpersample.loadSettingsFrom(settings);
    	m_bcfoutput.loadSettingsFrom(settings);
    	m_outbasepositions.loadSettingsFrom(settings);
    	m_outmapqual.loadSettingsFrom(settings);
    	m_strandbiaspval.loadSettingsFrom(settings);
    	m_uncompressedbcf.loadSettingsFrom(settings);
    	m_ifindex.loadSettingsFrom(settings);
    	m_usefaidxfile.loadSettingsFrom(settings);
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
    	m_ifreadsfile.validateSettings(settings);
    	m_excludereadsfile.validateSettings(settings);
    	m_capmapqual.validateSettings(settings);
    	m_ignorerg.validateSettings(settings);
    	m_minfrac.validateSettings(settings);
    	m_mingapreads.validateSettings(settings);
    	m_pileupregion.validateSettings(settings);
    	m_platformlist.validateSettings(settings);
    	
    	m_outpersample.validateSettings(settings);
    	m_bcfoutput.validateSettings(settings);
    	m_outbasepositions.validateSettings(settings);
    	m_outmapqual.validateSettings(settings);
    	m_strandbiaspval.validateSettings(settings);
    	m_uncompressedbcf.validateSettings(settings);
    	m_ifindex.validateSettings(settings);
    	m_usefaidxfile.validateSettings(settings);
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


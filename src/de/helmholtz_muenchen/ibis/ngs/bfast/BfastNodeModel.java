package de.helmholtz_muenchen.ibis.ngs.bfast;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;
import de.helmholtz_muenchen.ibis.utils.ngs.QSub;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;


public class BfastNodeModel extends NodeModel {

	//General parameters
	public static final String CFGKEY_THREADS = "threads";
	private final SettingsModelIntegerBounded m_threads = new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_THREADS, 4, 1, Integer.MAX_VALUE);
	
	
	//Ref. genome and indexing parameters
	public static final String CFGKEY_BFASTEXE = "bfast";
	public static final String CFGKEY_INPUTFASTA = "fasta";
	public static final String CFGKEY_MODELTYPE = "modeltype"; //if reg. genome is color or not
	public static final String CFGKEY_MASK = "mask"; //string of 1s and 0s with 1 at beginnning and end, e.g. 101
	public static final String CFGKEY_HASHWIDTH = "hashwidth"; //string of 1s and 0s with 1 at beginnning and end, e.g. 101
	
	public static final String CFGKEY_USESPLITDEPTH = "usesplitdepth";
	public static final String CFGKEY_SPLITDEPTH = "splitdepth";
	public static final String CFGKEY_USETMPDIR = "usetmpdir";
	public static final String CFGKEY_TMPDIR = "tmpdir";
	
	public static final String CFGKEY_USECONTIGS = "usecontigs";
	public static final String CFGKEY_USEPOS = "usepos";
	public static final String CFGKEY_USEEXONSFILE = "useexonsfile";
	public static final String CFGKEY_STARTCONTIG = "startcontig";
	public static final String CFGKEY_STARTPOS = "startpos";
	public static final String CFGKEY_ENDCONTIG = "endcontig";
	public static final String CFGKEY_ENDPOS = "endpos";
	public static final String CFGKEY_REPEATMASKER = "repeatmasker";
	public static final String CFGKEY_EXONSFILE = "exonsfile";
	
	
	private final SettingsModelString m_bfast = new SettingsModelString(
			BfastNodeModel.CFGKEY_BFASTEXE,"");
	private final SettingsModelString m_fasta = new SettingsModelString(
			BfastNodeModel.CFGKEY_INPUTFASTA,"");
	private final SettingsModelString m_modeltype = new SettingsModelString(
			BfastNodeModel.CFGKEY_MODELTYPE,"");
		private final SettingsModelString m_mask = new SettingsModelString(
			BfastNodeModel.CFGKEY_MASK,"1111111111111111111111 1111101110111010100101011011111 1011110101101001011000011010001111111 10111001101001100100111101010001011111 11111011011101111011111111 111111100101001000101111101110111 11110101110010100010101101010111111 111101101011011001100000101101001011101 1111011010001000110101100101100110100111 1111010010110110101110010110111011");
	private final SettingsModelIntegerBounded m_hashwidth = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_HASHWIDTH,12,1,Integer.MAX_VALUE);
	
	private final SettingsModelBoolean m_usesplitdepth = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USESPLITDEPTH, false);
	private final SettingsModelIntegerBounded m_splitdepth = new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_SPLITDEPTH, 0, 0, Integer.MAX_VALUE);
	private final SettingsModelBoolean m_usetmpdir = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USETMPDIR, false);
	private final SettingsModelString m_tmpdir = new SettingsModelString(
			BfastNodeModel.CFGKEY_TMPDIR,"");
	
	private final SettingsModelIntegerBounded m_startcontig = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_STARTCONTIG,1,1,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_endcontig = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_ENDCONTIG,1,1,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_startpos = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_STARTPOS,1,1,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_endpos = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_ENDPOS,1,1,Integer.MAX_VALUE);
	private final SettingsModelBoolean m_repeatmasker = new SettingsModelBoolean(BfastNodeModel.CFGKEY_REPEATMASKER, false);
	private final SettingsModelString m_exonsfile = new SettingsModelString(BfastNodeModel.CFGKEY_EXONSFILE,"");
	private final SettingsModelBoolean m_usecontigs = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USECONTIGS, false);
	private final SettingsModelBoolean m_usepos = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEPOS, false);
	private final SettingsModelBoolean m_useexonsfile = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEEXONSFILE, false);
	
	//Candidate Alignment Locations (CALs) parameters
	//public static final String CFGKEY_READSINPUTFILE = "readsinputfile";
	public static final String CFGKEY_COMPRESSTYPE = "compresstype";
	public static final String CFGKEY_LOADINDEXES = "loadindexes";
	public static final String CFGKEY_STRAND = "strand";
	public static final String CFGKEY_USEMAININDEXES = "usemainindexes";
	public static final String CFGKEY_MAININDEXES = "mainindexes";
	public static final String CFGKEY_USESECINDEXES = "usesecindexes";
	public static final String CFGKEY_SECINDEXES = "secindexes"; //secondary indexes
	public static final String CFGKEY_USEOFFSETS = "useoffsets";
	public static final String CFGKEY_OFFSETS = "offsets";
	public static final String CFGKEY_STARTREADNUM = "startreadnum";
	public static final String CFGKEY_ENDREADNUM = "endreadnum";
	public static final String CFGKEY_USEKEYSIZE = "usekeysize";
	public static final String CFGKEY_KEYSIZE = "keysize";
	public static final String CFGKEY_MAXKEYMATCHES = "maxkeymatches";
	public static final String CFGKEY_KEYMISSFRACTION = "keymissfraction";
	public static final String CFGKEY_MAXNUMMATCHES = "maxnummatches";
	public static final String CFGKEY_USEMAXREADS = "usemaxreads";
	public static final String CFGKEY_MAXREADS = "maxreads";
	
	//private final SettingsModelString m_readsinputfile = new SettingsModelString(RunBfastNodeModel.CFGKEY_READSINPUTFILE,"");
	private final SettingsModelString m_compresstype = new SettingsModelString(
			BfastNodeModel.CFGKEY_COMPRESSTYPE,"");
	private final SettingsModelString m_loadindexes = new SettingsModelString(
			BfastNodeModel.CFGKEY_LOADINDEXES,"");
	private final SettingsModelString m_strand = new SettingsModelString(
			BfastNodeModel.CFGKEY_STRAND,"");
	private final SettingsModelBoolean m_usemainindexes = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEMAININDEXES, false);
	private final SettingsModelBoolean m_usesecindexes = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USESECINDEXES, false);
	private final SettingsModelBoolean m_useoffsets = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEOFFSETS, false);
	private final SettingsModelString m_mainindexes = new SettingsModelString(
			BfastNodeModel.CFGKEY_MAININDEXES,"");
	private final SettingsModelString m_secindexes = new SettingsModelString(
			BfastNodeModel.CFGKEY_SECINDEXES,"");
	private final SettingsModelString m_offsets = new SettingsModelString(
			BfastNodeModel.CFGKEY_OFFSETS,"");
	private final SettingsModelIntegerBounded m_startreadnum = new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_STARTREADNUM, 1, 1, Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_endreadnum = new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_ENDREADNUM, 2147483647, 1, 2147483647);
	private final SettingsModelBoolean m_usekeysize = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEKEYSIZE, false);
	private final SettingsModelIntegerBounded m_keysize = new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_KEYSIZE,1,1,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_maxkeymatches = new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_MAXKEYMATCHES, 8, 1, Integer.MAX_VALUE);	
	private final SettingsModelDouble m_keymissfraction = new SettingsModelDouble(BfastNodeModel.CFGKEY_KEYMISSFRACTION, 1.0);
	private final SettingsModelIntegerBounded m_maxnummatches = new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_MAXNUMMATCHES, 384, 1, Integer.MAX_VALUE);
	private final SettingsModelBoolean m_usemaxreads = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEMAXREADS, false);
	private final SettingsModelIntegerBounded m_maxreads = new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_MAXREADS, 500000, 1, Integer.MAX_VALUE);	
	
	
	/**************************************
	*		Alignment&postprocessing	  *
	**************************************/
	//Align
	public static final String CFGKEY_GAPPEDALIGN = "gappedalign";
	public static final String CFGKEY_MASKCONSTRAINTS = "maskconstraints";
	public static final String CFGKEY_READSTART = "readstart";
	public static final String CFGKEY_READSTOP = "readstop";
	public static final String CFGKEY_OFFSET = "offset";
	public static final String CFGKEY_MAXMATCHES = "maxmatches";
	public static final String CFGKEY_AVGMISQUAL = "avgmisqual";
	
	//Postprocess
	public static final String CFGKEY_ALGO = "algo";
	public static final String CFGKEY_PAIRING = "pairing";
	public static final String CFGKEY_MINMAPQUAL = "minmapqual";
	public static final String CFGKEY_MINNORMSCORE = "minnormscore";
	public static final String CFGKEY_USEAVGDEV = "useavgdev";
	public static final String CFGKEY_INSSIZEAVG = "sizeavg";
	public static final String CFGKEY_INSSTDDEV = "stddev";

	//Align
	public static final int DEFAULT_READSTART = 1;
	public static final int DEFAULT_READSTOP = Integer.MAX_VALUE;
	public static final int DEFAULT_OFFSET = 20;
	public static final int DEFAULT_MAXMATCHES = 384;
	public static final int DEFAULT_AVGMISQUAL = 10;
	
	//Postprcess
	public static final String DEFAULT_ALGO = "no filtering";
	public static final String DEFAULT_PAIRING = "no pairing";
	public static final int DEFAULT_MINMAPQUAL = Integer.MIN_VALUE;
	public static final int DEFAULT_MINNORMSCORE = Integer.MIN_VALUE;
	public static final int DEFAULT_INSSIZEAVG = 0;
	public static final int DEFAULT_INSSTDDEV = 1;
	
	//Align
	private final SettingsModelString m_gappedalign = new SettingsModelString(
			BfastNodeModel.CFGKEY_GAPPEDALIGN,"");
	private final SettingsModelString m_maskconstraints = new SettingsModelString(
			BfastNodeModel.CFGKEY_MASKCONSTRAINTS,"");
	/*private final SettingsModelIntegerBounded m_readstart = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_READSTART,
			BfastNodeModel.DEFAULT_READSTART,
			0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_readstop = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_READSTOP,
			BfastNodeModel.DEFAULT_READSTOP,
			0,Integer.MAX_VALUE);*/
	private final SettingsModelIntegerBounded m_offset = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_OFFSET,
			BfastNodeModel.DEFAULT_OFFSET,
			0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_maxmatches = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_MAXMATCHES,
			BfastNodeModel.DEFAULT_MAXMATCHES,
			0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_avgmisqual = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_AVGMISQUAL,
			BfastNodeModel.DEFAULT_AVGMISQUAL,
			0,Integer.MAX_VALUE);
	//Postprocess
	private final SettingsModelString m_algo = new SettingsModelString(
			BfastNodeModel.CFGKEY_ALGO,DEFAULT_ALGO);
	private final SettingsModelString m_pairing = new SettingsModelString(
			BfastNodeModel.CFGKEY_PAIRING,DEFAULT_PAIRING);
	private final SettingsModelIntegerBounded m_minmapqual = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_MINMAPQUAL,
			BfastNodeModel.DEFAULT_MINMAPQUAL,
			Integer.MIN_VALUE,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_minnormscore = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_MINNORMSCORE,
			BfastNodeModel.DEFAULT_MINNORMSCORE,
			Integer.MIN_VALUE,Integer.MAX_VALUE);
	
	private final SettingsModelBoolean m_useavgdev = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEAVGDEV, false);
	private final SettingsModelIntegerBounded m_inssizeavg = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_INSSIZEAVG,
			BfastNodeModel.DEFAULT_INSSIZEAVG,
			0,Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_insstddev = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_INSSTDDEV,
			BfastNodeModel.DEFAULT_INSSTDDEV,
			0,Integer.MAX_VALUE);
	
	/**
     * Constructor for the node model.
     */
    protected BfastNodeModel() {
    
        super(1, 1);
        m_keysize.setEnabled(false);
        m_inssizeavg.setEnabled(false);
        m_insstddev.setEnabled(false);
        m_mainindexes.setEnabled(false);
        m_secindexes.setEnabled(false);
        m_offsets.setEnabled(false);
        m_splitdepth.setEnabled(false);
        m_tmpdir.setEnabled(false);
        m_startcontig.setEnabled(false);
        m_endcontig.setEnabled(false);
        m_startpos.setEnabled(false);
        m_endpos.setEnabled(false);
        m_exonsfile.setEnabled(false);
        m_maxreads.setEnabled(false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	

    	/**Initialize logfile**/
    	String folder = m_fasta.getStringValue();
    	folder = folder.substring(0,folder.lastIndexOf("/")+1);
    	String logfile = folder +"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("Bfast"));
    	/**logfile initialized**/
    	 
    	
    	 String path2readFile = inData[0].iterator().next().getCell(0).toString();
    	 String threads = " -n "+m_threads.getIntValue()+" ";
    	 String outBaseName = path2readFile.substring(path2readFile.lastIndexOf("/")+1,path2readFile.lastIndexOf("."));
    	 
    	 //Create Reference Genome   
    	 String model = "0";
    	 //-A 0 => no color ref. genome, -A 1 => color ref. genome
    	 if (m_modeltype.getStringValue().equals("Color space")) {
    		model = "1";
         }
    		 
    	 String call_1 = m_bfast.getStringValue() + " fasta2brg -f " + m_fasta.getStringValue() + " -A " + model; //fasta2brg creates ref. genome
    	     	
    	// begin QueueSub #################################################
 		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
 			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_Bfast-fasta2brg";
 			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
 			String logfle = path2readFile.substring(0,path2readFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
 			new QSub(call_1, name, memory, logfle, true);
 			logBuffer.append("QSub: " + call_1 + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
 		// end QueueSub ###################################################
 		} else {
	    	 Process p = Runtime.getRuntime().exec(call_1);
	    	 p.waitFor();
	    	 logBuffer.append(ShowOutput.getLogEntry(p, call_1));
 		}
    	 
    
    	String tmpdir = "";
    	if(m_usetmpdir.getBooleanValue()){
    		tmpdir = " -T "+m_tmpdir.getStringValue()+"/ ";
    	}
    	String splitdepth = "";
    	if(m_usesplitdepth.getBooleanValue()){
    		splitdepth = " -d "+m_splitdepth.getIntValue()+" ";
    	}
    	
    	//Create Indexed Reference Genome
    	String maskrepeats = "";
    	String positions = "";
    	String contigs = "";
    	String exonsfile = "";
    	if(m_repeatmasker.getBooleanValue()){
    		maskrepeats = " -R ";
    	}
    	if(m_usecontigs.getBooleanValue()){
    		contigs = " -s "+m_startcontig.getIntValue()+" -e "+m_endcontig.getIntValue()+" ";
    	}
    	if(m_usepos.getBooleanValue()){
    		positions = " -S "+m_startpos.getIntValue()+" -E "+m_endpos.getIntValue()+" ";
    	}
    	if(m_useexonsfile.getBooleanValue()){
    		exonsfile = " -x "+m_exonsfile.getStringValue()+" ";
    	}
    	    	
    	String call_2 = m_bfast.getStringValue() + " index -f " + m_fasta.getStringValue() + " -m " + m_mask.getStringValue() + " -w " + m_hashwidth.getIntValue() + tmpdir + splitdepth + contigs + positions + exonsfile + maskrepeats + threads+" > /dev/null 2>&1";
    	// begin QueueSub #################################################
 		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
 			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_Bfast-index";
 			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
 			String logfle = path2readFile.substring(0,path2readFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
 			new QSub(call_2, name, memory, logfle, true);
 			logBuffer.append("QSub: " + call_2 + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
 		// end QueueSub ###################################################
 		} else {
	    	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", call_2);
	    	/*for(Integer i=0; i<b.command().size();i++){
	    	System.out.println(b.command().get(i));
	    	}*/
	    	Process z;
	    	z = b.start();
	    	z.waitFor();
	    	logBuffer.append(ShowOutput.getLogEntry(z, call_2));
 		}

    	 
    	//Find Candidate Alignments
    	String compression = "";
    	
    	if(m_compresstype.getStringValue().equals("bz2")){
    		compression = "--bz2";
    	}
    	else if(m_compresstype.getStringValue().equals("gz")){
    		compression = "--gz";
    	}
    	String loadIndexes = "";
    	if(m_loadindexes.getStringValue().equals("yes")){
    		loadIndexes = "-l";
    	}
    	
    	String strand = "";
    	if(m_strand.getStringValue().equals("both")){
    		strand = "-w 0";   		
    	}
    	else if(m_strand.getStringValue().equals("forward")){
    		strand = "-w 1";
    	}
    	else if(m_strand.getStringValue().equals("reverse")){
    		strand = "-w 2";
    	}
    	
    	String startreadnum = "-s " + m_startreadnum.getIntValue();
    	String endreadnum = "-e " + m_endreadnum.getIntValue();
    	String maxkeymatches = "-K " + m_maxkeymatches.getIntValue();
    	String keymissfraction = "-F " + m_keymissfraction.getDoubleValue();
    	String maxnummatches = "-M " + m_maxnummatches.getIntValue();
    	String mainindexes = "";
    	String secindexes = "";
    	String offsets = "";
    	String keysize = "";
    	String maxreads = "";
    	if(m_usekeysize.getBooleanValue()){
    		keysize = "-k " + m_keysize.getIntValue();
    	}
    	if(m_usemainindexes.getBooleanValue() && !m_mainindexes.getStringValue().equals("")){
    		mainindexes = "-i "+m_mainindexes.getStringValue();
    	}
    	if(m_usesecindexes.getBooleanValue() && !m_secindexes.getStringValue().equals("")){
    		secindexes = "-I "+m_secindexes.getStringValue();
    	}
    	if(m_useoffsets.getBooleanValue() && !m_offsets.getStringValue().equals("")){
    		offsets = "-o "+m_offsets.getStringValue();
    	}
    	if(m_usemaxreads.getBooleanValue()){
    		maxreads = " -Q "+m_maxreads.getIntValue()+" ";
    	}
    	
    	String call_3 = m_bfast.getStringValue() + " match -f "+m_fasta.getStringValue()+" -r "+path2readFile+" "+
    		loadIndexes+" "+compression+" "+strand+" "+startreadnum+" "+endreadnum+" "+maxkeymatches+" "+keymissfraction+" "+maxnummatches+" "+
    		mainindexes+" "+secindexes+" "+offsets+" "+keysize+ maxreads +" -t" + threads;
	 	
    	//make path for output
    	String outPath = m_fasta.getStringValue().substring(0,m_fasta.getStringValue().lastIndexOf('/')+1);
    	outPath += "CAL_result";
    	call_3 += " > "+outPath;
    	
    	// begin QueueSub #################################################
 		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
 			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_Bfast-match";
 			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
 			String logfle = path2readFile.substring(0,path2readFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
 			new QSub(call_3, name, memory, logfle, true);
 			logBuffer.append("QSub: " + call_3 + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
 		// end QueueSub ###################################################
 		} else {
	    	//System.out.println(call_3+" >"+outPath);
	    	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", call_3);
			Process r = b.start();
			r.waitFor();
			logBuffer.append(ShowOutput.getLogEntry(r, call_3));
 		}
        
		/**************************************
		*		Alignment&postprocessing	  *
		**************************************/
		
		String path2bfast = m_bfast.getStringValue();
    	String path2reffile = m_fasta.getStringValue();
    	String path2match = outPath;
    	String color = model;

    	
    	String basePath = path2reffile.substring(0,path2reffile.lastIndexOf('/')+1);
    	String path2Alignoutfile = basePath+outBaseName+"_localalign.out";
    	
    	String gap="";
    	if(m_gappedalign.getStringValue().equals("ungapped")){
    		gap="-u ";
    	}
    	String mask="";
    	if(m_mask.getStringValue().equals("No")){
    		mask=" -U ";
    	}
    	String col=" -A "+color;
    	//String start = " -s "+m_readstart.getIntValue();
    	//String stop = " -e "+m_readstop.getIntValue();
    	String offset = " -o "+m_offset.getIntValue();
    	String maxmatch = " -M "+m_maxmatches.getIntValue();
    	String avgmisqual = " -q "+m_avgmisqual.getIntValue();
    	
    	String call = path2bfast+" localalign -f "+path2reffile+" -m "+path2match+" "+gap+mask+col+offset+maxmatch+avgmisqual + maxreads + threads +"> "+path2Alignoutfile;
    	//System.out.println(call);
    	// begin QueueSub #################################################
 		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
 			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_Bfast-localalign";
 			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
 			String logfle = path2readFile.substring(0,path2readFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
 			new QSub(call, name, memory, logfle, true);
 			logBuffer.append("QSub: " + call + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
 		// end QueueSub ###################################################
 		} else {
	    	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", call);
	    	Process p1 = b.start();
	    	p1.waitFor();
	        logBuffer.append(ShowOutput.getLogEntry(p1, call));
 		}
     
        //Begin Postprocessing
        int algo2use=-1;
        if(m_algo.getStringValue().equals("No filtering")){
        	algo2use=0;
        }else if(m_algo.getStringValue().equals("All alignments")){
        	algo2use=1;
        }else if(m_algo.getStringValue().equals("Uniquely aligned reads only")){
        	algo2use=2;
        }else if(m_algo.getStringValue().equals("Unique alignment with best score")){
        	algo2use=3;
        }else if(m_algo.getStringValue().equals("All alignments with best score")){
        	algo2use=4;
        }
        
        //System.out.print(m_pairing.getStringValue());
        int pairing2use=-1;
        if(m_pairing.getStringValue().equals("paired ends")){
        	pairing2use=0;
        }else if(m_pairing.getStringValue().equals("mate pairs")){
        	pairing2use=1;
        }else if(m_pairing.getStringValue().equals("no pairing")){
        	pairing2use=2;
         }
        
        String algo = " -a "+algo2use;
        String pairing = " -Y "+pairing2use;
        String minqual = " -m "+m_minmapqual.getIntValue();
        String minnorm = " -M "+m_minnormscore.getIntValue();
        String inssize = "";
        String insstd = "";
        if(m_useavgdev.getBooleanValue()){
        	inssize = " -v "+m_inssizeavg.getIntValue();
        	insstd = " -s "+m_insstddev.getIntValue();
        }
    	String path2Postoutfile = basePath+outBaseName+"_post.sam";
        
        String call1 = path2bfast+" postprocess -f "+path2reffile+" -i "+path2Alignoutfile+algo+pairing+minqual+minnorm+inssize+insstd + maxreads + threads + "> "+path2Postoutfile;
        // begin QueueSub #################################################
 		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
 			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_Bfast-postprocess";
 			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
 			String logfle = path2readFile.substring(0,path2readFile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
 			new QSub(call1, name, memory, logfle, true);
 			logBuffer.append("QSub: " + call1 + "\n");
			logBuffer.append("See external logfile: " + logfle + "\n");
 		// end QueueSub ###################################################
 		} else {
	        //System.out.println(call1);
	    	ProcessBuilder b1 = new ProcessBuilder("/bin/sh", "-c", call1);
	    	Process p2 = b1.start();
	    	p2.waitFor();
	        logBuffer.append(ShowOutput.getLogEntry(p2, call1));
 		}
        
        
        
        //Create Output
        
        DataColumnSpecCreator col1 = new DataColumnSpecCreator("Path2SAMFile", StringCell.TYPE);
        DataColumnSpecCreator col2 = new DataColumnSpecCreator("Path2RefFile", StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec(),col2.createSpec()};
     
    	DataTableSpec table = new DataTableSpec(cols);
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	StringCell cl1 = new StringCell(path2Postoutfile);
    	StringCell cl2 = new StringCell(path2reffile);
    	DataCell[] c = new DataCell[]{cl1,cl2};
    	

    	DefaultRow d = new DefaultRow("Row0",c);
    	cont.addRowToTable(d);
    	cont.close();
    	BufferedDataTable out = cont.getTable();

    	pushFlowVariableString("BAMSAMINFILE",path2Postoutfile);
    	
    	logBuffer.append(ShowOutput.getNodeEndTime());
    	ShowOutput.writeLogFile(logBuffer);
    	
        return new BufferedDataTable[]{out};

    }

    /**
     * {@inheritDoc}
     */
    protected void reset() {
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {

     	if(m_startreadnum.getIntValue() > m_endreadnum.getIntValue()){
            throw new InvalidSettingsException("Start read number can't be bigger than end read number!");
    	}

     	if(m_fasta.getStringValue().length() > 1) {
	     	if(!FileValidator.checkFastaFormat(m_fasta.getStringValue())){
	     		throw new InvalidSettingsException("Given file is not in fasta format!");
	     	}
     	}
     	
     	
     	// Check if single-end
    	if(getAvailableInputFlowVariables().get("readType").getStringValue().equals("paired-end")) {
    		setWarningMessage("Paired-end mapping is not implemented yet. Bfast will use single-end mapping.");
    	}
    	
    	// Check if FastQ file
    	String isBam = getAvailableInputFlowVariables().get("isBAM").getStringValue();
    	if(isBam.equals("true")) {
    			throw new InvalidSettingsException("Bfast does not support BAM files. Please choose a FastQ or FastA file containing your reads.");
    	}
     	
     	//check if masks have "1" at start and end
     	String mask = m_mask.getStringValue();
     	String[] masks = mask.split(" ");
     	for(int i=0;i<masks.length;i++){
     		String curr_mask = masks[i];
     		char first = curr_mask.charAt(0);
         	char last = curr_mask.charAt(curr_mask.length()-1);
         	if( !(first=='1' && last=='1') ){
         		throw new InvalidSettingsException("Invalid mask. Mask has to start and end with '1'");	
         	
         	}

     	}
        //Version control
     	if(m_bfast.getStringValue().length()>0){
            if(FileValidator.versionControl(m_bfast.getStringValue(),"BFAST")==1){
            	setWarningMessage("WARNING: You are using a newer BFAST version than "+FileValidator.BFAST_VERSION +"! This may cause problems");
            }else if(FileValidator.versionControl(m_bfast.getStringValue(),"BFAST")==2){
            	throw new InvalidSettingsException("You are using a outdated version of BFAST! Please update your version");
            }else if(FileValidator.versionControl(m_bfast.getStringValue(),"BFAST")==-1){
            	System.out.println("Something wrong here");
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
    	m_bfast.saveSettingsTo(settings);
    	m_fasta.saveSettingsTo(settings);
    	m_modeltype.saveSettingsTo(settings);
    	m_hashwidth.saveSettingsTo(settings);
    	m_mask.saveSettingsTo(settings);
    	m_threads.saveSettingsTo(settings);
 	
    	//m_readsinputfile.saveSettingsTo(settings);
    	m_compresstype.saveSettingsTo(settings);
    	m_loadindexes.saveSettingsTo(settings);
    	m_strand.saveSettingsTo(settings);
    	
    	m_mainindexes.saveSettingsTo(settings);
    	m_secindexes.saveSettingsTo(settings); //secondary indexes
    	m_offsets.saveSettingsTo(settings);
    	m_startreadnum.saveSettingsTo(settings);
    	m_endreadnum.saveSettingsTo(settings);
    	m_keysize.saveSettingsTo(settings);
    	m_maxkeymatches.saveSettingsTo(settings);
    	m_keymissfraction.saveSettingsTo(settings);
    	m_maxnummatches.saveSettingsTo(settings);
        m_usemainindexes.saveSettingsTo(settings);
        m_usesecindexes.saveSettingsTo(settings);
        m_useoffsets.saveSettingsTo(settings);
    	
    	m_usesplitdepth.saveSettingsTo(settings);;
    	m_splitdepth.saveSettingsTo(settings);;
    	m_usetmpdir.saveSettingsTo(settings);
    	m_tmpdir.saveSettingsTo(settings);
    	
    	m_startcontig.saveSettingsTo(settings);
    	m_endcontig.saveSettingsTo(settings);
    	m_startpos.saveSettingsTo(settings);
    	m_endpos.saveSettingsTo(settings);
    	m_repeatmasker.saveSettingsTo(settings);
    	m_exonsfile.saveSettingsTo(settings);

    	m_usecontigs.saveSettingsTo(settings);
    	m_usepos.saveSettingsTo(settings);
    	m_useexonsfile.saveSettingsTo(settings);
    	
    	m_usemaxreads.saveSettingsTo(settings);
    	m_maxreads.saveSettingsTo(settings);
    	
    	//Align&postprocess
        m_avgmisqual.saveSettingsTo(settings);
        m_gappedalign.saveSettingsTo(settings);
        m_maskconstraints.saveSettingsTo(settings);
        m_maxmatches.saveSettingsTo(settings);
        m_offset.saveSettingsTo(settings);
        //m_readstart.saveSettingsTo(settings);
        //m_readstop.saveSettingsTo(settings);
        m_algo.saveSettingsTo(settings);
        m_pairing.saveSettingsTo(settings);
        m_minmapqual.saveSettingsTo(settings);
        m_minnormscore.saveSettingsTo(settings);
        m_useavgdev.saveSettingsTo(settings);
        m_inssizeavg.saveSettingsTo(settings);
        m_insstddev.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_bfast.loadSettingsFrom(settings);
    	m_fasta.loadSettingsFrom(settings);
    	m_modeltype.loadSettingsFrom(settings);
    	m_hashwidth.loadSettingsFrom(settings);
    	m_mask.loadSettingsFrom(settings);
    	m_threads.loadSettingsFrom(settings);
    	
    	//m_readsinputfile.loadSettingsFrom(settings);
    	m_compresstype.loadSettingsFrom(settings);
    	m_loadindexes.loadSettingsFrom(settings);
    	m_strand.loadSettingsFrom(settings);
    	
    	m_mainindexes.loadSettingsFrom(settings);
    	m_secindexes.loadSettingsFrom(settings); //secondary indexes
    	m_offsets.loadSettingsFrom(settings);
    	m_startreadnum.loadSettingsFrom(settings);
    	m_endreadnum.loadSettingsFrom(settings);
    	m_keysize.loadSettingsFrom(settings);
    	m_maxkeymatches.loadSettingsFrom(settings);
    	m_keymissfraction.loadSettingsFrom(settings);
    	m_maxnummatches.loadSettingsFrom(settings);
        m_usemainindexes.loadSettingsFrom(settings);
        m_usesecindexes.loadSettingsFrom(settings);
        m_useoffsets.loadSettingsFrom(settings);
        
    	m_usesplitdepth.loadSettingsFrom(settings);
    	m_splitdepth.loadSettingsFrom(settings);
    	m_usetmpdir.loadSettingsFrom(settings);
    	m_tmpdir.loadSettingsFrom(settings);
    	m_startcontig.loadSettingsFrom(settings);
    	m_endcontig.loadSettingsFrom(settings);
    	m_startpos.loadSettingsFrom(settings);
    	m_endpos.loadSettingsFrom(settings);
    	m_repeatmasker.loadSettingsFrom(settings);
    	m_exonsfile.loadSettingsFrom(settings);
    	m_usecontigs.loadSettingsFrom(settings);
    	m_usepos.loadSettingsFrom(settings);
    	m_useexonsfile.loadSettingsFrom(settings);
    	m_usemaxreads.loadSettingsFrom(settings);
    	m_maxreads.loadSettingsFrom(settings);
        //Align&postprocess
        m_avgmisqual.loadSettingsFrom(settings);
        m_gappedalign.loadSettingsFrom(settings);
        m_maskconstraints.loadSettingsFrom(settings);
        m_maxmatches.loadSettingsFrom(settings);
        m_offset.loadSettingsFrom(settings);
        //m_readstart.loadSettingsFrom(settings);
        //m_readstop.loadSettingsFrom(settings);
        m_algo.loadSettingsFrom(settings);
        m_pairing.loadSettingsFrom(settings);
        m_minmapqual.loadSettingsFrom(settings);
        m_minnormscore.loadSettingsFrom(settings);
        m_inssizeavg.loadSettingsFrom(settings);
        m_insstddev.loadSettingsFrom(settings);
        m_useavgdev.loadSettingsFrom(settings);
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_bfast.validateSettings(settings);
    	m_fasta.validateSettings(settings);
    	m_modeltype.validateSettings(settings);
    	m_hashwidth.validateSettings(settings);
    	m_mask.validateSettings(settings);
    	m_threads.validateSettings(settings);
    	
    	//m_readsinputfile.validateSettings(settings);
    	m_compresstype.validateSettings(settings);
    	m_loadindexes.validateSettings(settings);
    	m_strand.validateSettings(settings);
    	
    	m_offsets.validateSettings(settings);
    	m_startreadnum.validateSettings(settings);
    	m_endreadnum.validateSettings(settings);
    	m_keysize.validateSettings(settings);
    	m_maxkeymatches.validateSettings(settings);
    	m_keymissfraction.validateSettings(settings);
    	m_maxnummatches.validateSettings(settings);
        m_usemainindexes.validateSettings(settings);
        m_usesecindexes.validateSettings(settings);
        m_useoffsets.validateSettings(settings);
        
    	m_usesplitdepth.validateSettings(settings);
    	m_splitdepth.validateSettings(settings);
    	m_usetmpdir.validateSettings(settings);
    	m_tmpdir.validateSettings(settings);
        
    	m_startcontig.validateSettings(settings);
    	m_endcontig.validateSettings(settings);
    	m_startpos.validateSettings(settings);
    	m_endpos.validateSettings(settings);
    	m_repeatmasker.validateSettings(settings);
    	m_exonsfile.validateSettings(settings);
    	m_usecontigs.validateSettings(settings);
    	m_usepos.validateSettings(settings);
    	m_useexonsfile.validateSettings(settings);
    	m_usemaxreads.validateSettings(settings);
    	m_maxreads.validateSettings(settings);
        //Align&postprocess
        m_avgmisqual.validateSettings(settings);
        m_gappedalign.validateSettings(settings);
        m_maskconstraints.validateSettings(settings);
        m_maxmatches.validateSettings(settings);
        m_offset.validateSettings(settings);
       // m_readstart.validateSettings(settings);
       // m_readstop.validateSettings(settings);
        m_algo.validateSettings(settings);
        m_pairing.validateSettings(settings);
        m_minmapqual.validateSettings(settings);
        m_minnormscore.validateSettings(settings);
        m_inssizeavg.validateSettings(settings);
        m_insstddev.validateSettings(settings);
        m_useavgdev.validateSettings(settings);

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
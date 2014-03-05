package de.helmholtz_muenchen.ibis.ngs.bfast;

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
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;


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
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(BfastNodeModel.class);
		
	//The Output Col Names
	public static final String OUT_COL1 = "Path2SAMFile";
	public static final String OUT_COL2 = "Path2RefFile";
	
	
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
    	 
    	ArrayList<String> command = new ArrayList<String>();
    	
    	
    	 String path2readFile = inData[0].iterator().next().getCell(0).toString();
    	 String threads = "-n "+m_threads.getIntValue();
    	 String outBaseName = path2readFile.substring(path2readFile.lastIndexOf("/")+1,path2readFile.lastIndexOf("."));
    	 
    	 //Create Reference Genome   
    	 String model = "0";
    	 //-A 0 => no color ref. genome, -A 1 => color ref. genome
    	 if (m_modeltype.getStringValue().equals("Color space")) {
    		model = "1";
         }
    	
    	//fasta2brg creates ref. genome
    	 command.add(m_bfast.getStringValue()+" fasta2brg");
    	 command.add("-f "+m_fasta.getStringValue());
    	 command.add("-A "+model);

     	/**
     	 * Execute
     	 */
     	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
     	logBuffer.append(ShowOutput.getNodeEndTime());
     	ShowOutput.writeLogFile(logBuffer);
     	command = new ArrayList<String>();	//Clear Array
    

    	
    	//Create Indexed Reference Genome
    	command.add(m_bfast.getStringValue()+" index");
    	
    	command.add("-f "+m_fasta.getStringValue());
    	command.add("-m "+m_mask.getStringValue());
    	command.add("-w " + m_hashwidth.getIntValue());
    	    	
    	if(m_usetmpdir.getBooleanValue()){
    		command.add("-T "+m_tmpdir.getStringValue()+"/");
    	}
    	
    	if(m_usesplitdepth.getBooleanValue()){
    		command.add("-d "+m_splitdepth.getIntValue());
    	}
    	
    	if(m_repeatmasker.getBooleanValue()){
    		command.add("-R ");
    	}
    	if(m_usecontigs.getBooleanValue()){
    		command.add("-s "+m_startcontig.getIntValue()+" -e "+m_endcontig.getIntValue());
    	}
    	if(m_usepos.getBooleanValue()){
    		command.add("-S "+m_startpos.getIntValue()+" -E "+m_endpos.getIntValue());
    	}
    	if(m_useexonsfile.getBooleanValue()){
    		command.add("-x "+m_exonsfile.getStringValue());
    	}
    	
    	command.add(threads);

     	/**
     	 * Execute
     	 */
     	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
     	logBuffer.append(ShowOutput.getNodeEndTime());
     	ShowOutput.writeLogFile(logBuffer);
     	command = new ArrayList<String>();	//Clear Array
    	 
     	
     	
//    	String call_3 = m_bfast.getStringValue() + " match -f "+m_fasta.getStringValue()+" -r "+path2readFile+" "+
//        		loadIndexes+" "+compression+" "+strand+" "+startreadnum+" "+endreadnum+" "+maxkeymatches+" "+keymissfraction+" "+maxnummatches+" "+
//        		mainindexes+" "+secindexes+" "+offsets+" "+keysize+ maxreads +" -t" + threads;
//    	 	
//        	//make path for output
//        	String outPath = m_fasta.getStringValue().substring(0,m_fasta.getStringValue().lastIndexOf('/')+1);
//        	outPath += "CAL_result";
//        	call_3 += " > "+outPath;

     	
    	command.add(m_bfast.getStringValue()+" match");    	
     	command.add("-f "+m_fasta.getStringValue());
     	command.add("-r "+path2readFile);
    	if(m_loadindexes.getStringValue().equals("yes")){
    		command.add("-l");
    	}
    	if(m_compresstype.getStringValue().equals("bz2")){
    		command.add("--bz2");
    	}
    	else if(m_compresstype.getStringValue().equals("gz")){
    		command.add("--gz");
    	}

    	if(m_strand.getStringValue().equals("both")){
    		command.add("-w 0");   		
    	}
    	else if(m_strand.getStringValue().equals("forward")){
    		command.add("-w 1");
    	}
    	else if(m_strand.getStringValue().equals("reverse")){
    		command.add("-w 2");
    	}
    	
    	command.add("-s " + m_startreadnum.getIntValue());
    	command.add("-e " + m_endreadnum.getIntValue());
    	command.add("-K " + m_maxkeymatches.getIntValue());
    	command.add("-F " + m_keymissfraction.getDoubleValue());
    	command.add("-M " + m_maxnummatches.getIntValue());

    	if(m_usekeysize.getBooleanValue()){
    		command.add("-k " + m_keysize.getIntValue());
    	}
    	if(m_usemainindexes.getBooleanValue() && !m_mainindexes.getStringValue().equals("")){
    		command.add("-i "+m_mainindexes.getStringValue());
    	}
    	if(m_usesecindexes.getBooleanValue() && !m_secindexes.getStringValue().equals("")){
    		command.add("-I "+m_secindexes.getStringValue());
    	}
    	if(m_useoffsets.getBooleanValue() && !m_offsets.getStringValue().equals("")){
    		command.add("-o "+m_offsets.getStringValue());
    	}
    	if(m_usemaxreads.getBooleanValue()){
    		command.add("-Q "+m_maxreads.getIntValue());
    	}
    	//command.add("-t" + threads);
	 	
    	//make path for output
    	String outPath = m_fasta.getStringValue().substring(0,m_fasta.getStringValue().lastIndexOf('/')+1);
    	outPath += "CAL_result";

     	/**
     	 * Execute
     	 */
     	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER,outPath);
     	logBuffer.append(ShowOutput.getNodeEndTime());
     	ShowOutput.writeLogFile(logBuffer);
     	command = new ArrayList<String>();	//Clear Array
     	

        
		/**************************************
		*		Alignment&postprocessing	  *
		**************************************/
		
		String path2bfast = m_bfast.getStringValue();
    	String path2reffile = m_fasta.getStringValue();
    	String path2match = outPath;
    	String color = model;

    	
    	String basePath = path2reffile.substring(0,path2reffile.lastIndexOf('/')+1);
    	String path2Alignoutfile = basePath+outBaseName+"_localalign.out";
    	


    	command.add(path2bfast+" localalign");
    	command.add("-f "+path2reffile);
    	command.add("-m "+path2match);
    	if(m_gappedalign.getStringValue().equals("ungapped")){
    		command.add("-u");
    	}
    	if(m_mask.getStringValue().equals("No")){
    		command.add("-U");
    	}
    	command.add("-A "+color);
    	command.add("-o "+m_offset.getIntValue());
    	command.add("-M "+m_maxmatches.getIntValue());
    	command.add("-q "+m_avgmisqual.getIntValue());
    	if(m_usemaxreads.getBooleanValue()){
    		command.add("-Q "+m_maxreads.getIntValue());
    	}
    	command.add(threads);
    	
     	/**
     	 * Execute
     	 */
     	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER,path2Alignoutfile);
     	logBuffer.append(ShowOutput.getNodeEndTime());
     	ShowOutput.writeLogFile(logBuffer);
     	command = new ArrayList<String>();	//Clear Array
    	
    	
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
        
    	String path2Postoutfile = basePath+outBaseName+"_post.sam";
        
    	command.add(path2bfast+" postprocess");
    	command.add("-f "+path2reffile);
    	command.add("-i "+path2Alignoutfile);
    	command.add("-a "+algo2use);
    	command.add("-Y "+pairing2use);
    	command.add("-m "+m_minmapqual.getIntValue());
    	command.add("-M "+m_minnormscore.getIntValue());

        if(m_useavgdev.getBooleanValue()){
        	command.add("-v "+m_inssizeavg.getIntValue());
        	command.add("-s "+m_insstddev.getIntValue());
        }
    	if(m_usemaxreads.getBooleanValue()){
    		command.add("-Q "+m_maxreads.getIntValue());
    	}
    	command.add(threads);
    	
     	/**
     	 * Execute
     	 */
     	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER,path2Postoutfile);
     	logBuffer.append(ShowOutput.getNodeEndTime());
     	ShowOutput.writeLogFile(logBuffer);
        
    	/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(path2Postoutfile),
    			(FileCell) FileCellFactory.create(path2reffile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	pushFlowVariableString("BAMSAMINFILE",path2Postoutfile);
        
    	
        return new BufferedDataTable[]{outTable};

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
package de.helmholtz_muenchen.ibis.ngs.samtools;

import java.io.File;
import java.io.FileWriter;
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
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;
import de.helmholtz_muenchen.ibis.utils.ngs.QSub;
import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;


/**
 * This is the model implementation of SamTools.
 * 
 */

public class SamToolsNodeModel extends NodeModel {
    
	
	protected static boolean useSamtools = true;
	protected static boolean useBamfile = true;
	protected static boolean useFastafile = true;
	
	public static final String CFGKEY_UTILITY = "utility";
	
	public static final String CFGKEY_SAMTOOLS = "samtools";
	public static final String CFGKEY_BAMFILE = "bamfile";
	public static final String CFGKEY_REFSEQFILE = "refseqfile";
	
	
	private final SettingsModelString m_utility = new SettingsModelString(SamToolsNodeModel.CFGKEY_UTILITY, "");
	private final SettingsModelString m_samtools = new SettingsModelString(SamToolsNodeModel.CFGKEY_SAMTOOLS, "");
	private final SettingsModelString m_bamfile = new SettingsModelString(SamToolsNodeModel.CFGKEY_BAMFILE, "");
	private final SettingsModelString m_refseqfile = new SettingsModelString(SamToolsNodeModel.CFGKEY_REFSEQFILE, "");
	

	/**Parameters for "fillmd"**/
	public static final String CFGKEY_CHANGEIDENTBASES = "changeIdentBases";
	public static final String CFGKEY_USECOMPRESSION = "useCompression";
	public static final String CFGKEY_COMPRESSION = "compression";
	public static final String CFGKEY_INPUTISSAM = "inputIsSam";
	public static final String CFGKEY_MODIFYQUAL = "modifyQual";
	public static final String CFGKEY_BQTAG = "bqTag";
	public static final String CFGKEY_EXTENDEDBAQ = "extendedBAQ";
	public static final String CFGKEY_DOCAPMAPQUAL = "doCapMapQual";
	public static final String CFGKEY_CAPMAPQUAL = "capMapQual";
	
	private final SettingsModelBoolean m_changeIdentBases = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_CHANGEIDENTBASES, false);
	private final SettingsModelBoolean m_useCompression = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USECOMPRESSION, false);	
	private final SettingsModelString m_compression = new SettingsModelString(SamToolsNodeModel.CFGKEY_COMPRESSION, "");
	private final SettingsModelBoolean m_inputIsSAM = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_INPUTISSAM, false);
	private final SettingsModelBoolean m_modifyQual = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MODIFYQUAL, false);
	private final SettingsModelBoolean m_bqTag = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_BQTAG, false);
	private final SettingsModelBoolean m_extendedBAQ = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_EXTENDEDBAQ, false);
	//private final SettingsModelBoolean m_doCapMapQual = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_DOCAPMAPQUAL, false);
	//private final SettingsModelIntegerBounded m_capMapQual = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_CAPMAPQUAL, 1, 1, Integer.MAX_VALUE);
	/**Parameters for rmdup**/
	public static final String CFGKEY_REMOVEDUP = "removeDup";
	public static final String CFGKEY_TREATPE = "treatPE";
	private final SettingsModelBoolean m_removeDup = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_REMOVEDUP, false);
	private final SettingsModelBoolean m_treatPE = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_TREATPE, false);
	/**Parameters for cat**/
	public static final String CFGKEY_USEHEADERSAM = "useHeaderSAM";
	public static final String CFGKEY_HEADERSAM = "headerSAM";
	public static final String CFGKEY_INBAM1 = "inBAM1";
	//public static final String CFGKEY_INBAM2 = "inBAM2";
	private final SettingsModelBoolean m_useHeaderSAM = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEHEADERSAM, false);
	private final SettingsModelString m_headerSAM = new SettingsModelString(SamToolsNodeModel.CFGKEY_HEADERSAM, "");
	private final SettingsModelString m_inBAM1 = new SettingsModelString(SamToolsNodeModel.CFGKEY_INBAM1, "");
	//private final SettingsModelString m_inBAM2 = new SettingsModelString(SamToolsNodeModel.CFGKEY_INBAM2, "");
	/**Parameters for reheader**/
	public static final String CFGKEY_REHINSAM = "rehInSAM";
	//public static final String CFGKEY_REHINBAM = "rehInBAM";
	//private final SettingsModelString m_rehInBAM = new SettingsModelString(SamToolsNodeModel.CFGKEY_REHINBAM, "");
	private final SettingsModelString m_rehInSAM = new SettingsModelString(SamToolsNodeModel.CFGKEY_REHINSAM, "");
	/**Parameters for merge**/
	public static final String CFGKEY_MINBAM1 = "minbam1";
	//public static final String CFGKEY_MINBAM2 = "minbam2";
	public static final String CFGKEY_MCOMPRESSION = "mcompression"; //zlib compression
	public static final String CFGKEY_MFORCE = "mforce";
	public static final String CFGKEY_USEMHFILE = "usemhfile";
	public static final String CFGKEY_MHFILE = "mhfile";
	public static final String CFGKEY_MSORTED = "msorted";
	public static final String CFGKEY_MREGION = "mregion";
	public static final String CFGKEY_USEMREGION = "usemregion";
	public static final String CFGKEY_MRGTAG = "mrgtag";
	public static final String CFGKEY_MUNCOMPRESSED = "muncompressed";
	private final SettingsModelString m_minbam1 = new SettingsModelString(SamToolsNodeModel.CFGKEY_MINBAM1, "");
	//private final SettingsModelString m_minbam2 = new SettingsModelString(SamToolsNodeModel.CFGKEY_MINBAM2, "");
	private final SettingsModelBoolean m_mcompression = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MCOMPRESSION, false);
	private final SettingsModelBoolean m_mforce = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MFORCE, false);
	private final SettingsModelBoolean m_usemhfile = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEMHFILE, false);
	private final SettingsModelString m_mhfile = new SettingsModelString(SamToolsNodeModel.CFGKEY_MHFILE, "");
	private final SettingsModelBoolean m_msorted = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MSORTED, false);
	private final SettingsModelString m_mregion = new SettingsModelString(SamToolsNodeModel.CFGKEY_MREGION, "");
	private final SettingsModelBoolean m_usemregion = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEMREGION, false);
	private final SettingsModelBoolean m_mrgtag = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MRGTAG, false);
	private final SettingsModelBoolean m_muncompressed = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MUNCOMPRESSED, false);
	/**faidx**/
	//public static final String CFGKEY_INFASTA = "infasta";
	//private final SettingsModelString m_infasta = new SettingsModelString(SamToolsNodeModel.CFGKEY_INFASTA, "");
	/**phase**/
	public static final String CFGKEY_BLOCKLENGTH = "blocklength";
	public static final String CFGKEY_PREFIX = "prefix";
	public static final String CFGKEY_HETPHRED = "hetphred";
	public static final String CFGKEY_MINQUAL = "minqual";
	public static final String CFGKEY_MAXDEPTH = "maxdepth";
	public static final String CFGKEY_FIXCHIMERAS = "fixchimeras";
	//public static final String CFGKEY_DROPAMBIG = "dropambig"; //does not work (via console) for some reason, maybe deprecated
	private final SettingsModelIntegerBounded m_blocklength = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_BLOCKLENGTH, 13, 1, Integer.MAX_VALUE);
	private final SettingsModelString m_prefix = new SettingsModelString(SamToolsNodeModel.CFGKEY_PREFIX,"phase");
	private final SettingsModelIntegerBounded m_hetphred = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_HETPHRED, 37, 1, Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_minqual = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_MINQUAL, 13, 1, Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_maxdepth = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_MAXDEPTH, 256, 1, Integer.MAX_VALUE);
	private final SettingsModelBoolean m_fixchimeras = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_FIXCHIMERAS, false);
	//private final SettingsModelBoolean m_dropambig = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_DROPAMBIG, false);
	
	protected static Boolean optionalPort = false;
	
	/**
     * Constructor for the node model.
     */
    protected SamToolsNodeModel() {
    
        super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(0));
       
        //m_refseqfile.setEnabled(false);
        //fillmd
        m_changeIdentBases.setEnabled(false);
    	m_compression.setEnabled(false);
    	m_useCompression.setEnabled(false);
    	m_inputIsSAM.setEnabled(false);
    	m_modifyQual.setEnabled(false);
    	m_bqTag.setEnabled(false);
    	m_extendedBAQ.setEnabled(false);
    	//m_doCapMapQual.setEnabled(false);
    	//m_capMapQual.setEnabled(false);
    	//mdup
    	m_removeDup.setEnabled(false);
    	m_treatPE.setEnabled(false);
    	//cat
    	m_headerSAM.setEnabled(false);
    	//m_inBAM1.setEnabled(false);
    	//m_inBAM2.setEnabled(false);
    	//reheader
    	m_rehInSAM.setEnabled(false);
    	//m_rehInBAM.setEnabled(false);
    	//merge
    	m_mcompression.setEnabled(false);
    	m_mforce.setEnabled(false);
    	m_usemhfile.setEnabled(false);
    	m_mhfile.setEnabled(false);
    	m_msorted.setEnabled(false);
    	m_mregion.setEnabled(false);
    	m_usemregion.setEnabled(false);
    	m_mrgtag.setEnabled(false);
    	m_muncompressed.setEnabled(false);
    	//m_minbam1.setEnabled(false);
    	//m_minbam2.setEnabled(false);
    	//faidx
    	//m_infasta.setEnabled(false);
    	//phase
    	m_blocklength.setEnabled(false);
    	m_prefix.setEnabled(false);
    	m_hetphred.setEnabled(false);
    	m_minqual.setEnabled(false);
    	m_maxdepth.setEnabled(false);
    	m_fixchimeras.setEnabled(false);
    	//m_dropambig.setEnabled(false);
    }

    /**
     * {@inheritDoc}
     * @return 
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	

    	
    	String path2samtools = "";
    	String path2bamfile = "";
    	String path2seqfile = "";
    	
    	if(optionalPort){
    		path2samtools = inData[0].iterator().next().getCell(0).toString();
    		path2bamfile = inData[0].iterator().next().getCell(1).toString();
    		//path2seqfile = inData[0].iterator().next().getCell(2).toString();
    		path2seqfile = m_refseqfile.getStringValue();
    	}
    	else{
    		path2samtools = m_samtools.getStringValue();
    		path2bamfile = m_bamfile.getStringValue();
    		path2seqfile = m_refseqfile.getStringValue();
    	}
    	
    	/**Initialize logfile**/
    	String folder = "";
    	if(path2seqfile.length()>0){
    		folder = path2seqfile;
    	}
    	else{
    		folder = path2bamfile;
    	}
    	folder = folder.substring(0,folder.lastIndexOf("/")+1);
    	String logfile = folder +"logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("SamTools"));
    	/**logfile initialized**/
    	
    	String basePathWithFileName = "";
    	if(path2bamfile != ""){
    	   	basePathWithFileName = path2bamfile.substring(0,path2bamfile.lastIndexOf("."));	
    	}
    	else{
    		basePathWithFileName = path2seqfile.substring(0,path2seqfile.lastIndexOf("."));	
    	}
 
    	String utility = m_utility.getStringValue();
    	String com = "";
    	

    	
    	if(utility.equals("flagstat")) {
    		com = path2samtools + " " + utility + " " + path2bamfile + " > " + basePathWithFileName + "_" + utility + ".txt";
    	} else if(utility.equals("idxstats")) {
    		String outFile = basePathWithFileName + "_" + utility + ".txt";
    		com = path2samtools + " " + utility + " " + path2bamfile + " >> " + outFile;
    		File file = new File(outFile);
    		FileWriter writer = new FileWriter(file ,false);
    		writer.write("sequence name\tsequence length\t# mapped reads\t# unmapped reads\n");
    		writer.flush();
    		writer.close();
    	} else if(utility.equals("depth")) {
    		com = path2samtools + " " + utility + " " + path2bamfile + " > " + basePathWithFileName + "_" + utility + ".txt";

    	} 
    	
    	else if(utility.equals("fillmd")) {
    		//com = path2samtools + " " + utility + " " + path2bamfile + " > " + basePathWithFileName + "_" + utility + ".txt";
    		String outFileExtension = ".sam";
    		StringBuffer arguments = new StringBuffer(5);
    		if(m_changeIdentBases.getBooleanValue()){
    			arguments.append("-e ");
    		}
    		if(m_useCompression.getBooleanValue()){
    			outFileExtension = ".bam";
    			if(m_compression.getStringValue().equals("uncompressed")){
    				arguments.append("-u ");
    			}
    			else{
    				arguments.append("-b ");
    			}
    		}
    		if(m_inputIsSAM.getBooleanValue()){
    			arguments.append("-S ");
    		}
    		if(m_modifyQual.getBooleanValue()){
    			arguments.append("-A ");
    		}
    		if(m_bqTag.getBooleanValue()){
    			arguments.append("-r ");
    		}
    		if(m_extendedBAQ.getBooleanValue()){
    			arguments.append("-E ");
    		}
    		/*if(m_doCapMapQual.getBooleanValue()){
    			String value = " "+m_capMapQual.getIntValue();
    			arguments.append("-I" + value+" ");
    		}*/
    		
    		com = path2samtools + " " + utility + " " + arguments.toString() + path2bamfile + " " + path2seqfile + " > " + basePathWithFileName + "_" + utility + outFileExtension;
    	} 
    	else if (utility.equals("fixmate")) {
    		com = path2samtools + " " + utility + " "  + path2bamfile + " " + basePathWithFileName + "_" + utility + ".bam";
    	}
    	else if (utility.equals("rmdup")) {
    		String params = "";
    		if(m_removeDup.getBooleanValue()){
    			params += "-s ";
    		}
    		if(m_treatPE.getBooleanValue()){
    			params += "-S ";
    		}
    		com = path2samtools + " " + utility + " "  + params + path2bamfile + " " + basePathWithFileName + "_" + utility + ".bam";
    	}
    	else if (utility.equals("cat")) {
    		String outfile = path2bamfile.substring(0,path2bamfile.lastIndexOf("."));
    		outfile += "_" + m_inBAM1.getStringValue().substring(m_inBAM1.getStringValue().lastIndexOf("/")+1,m_inBAM1.getStringValue().lastIndexOf("."));
    		outfile += "_" + utility + ".bam";
    		String params = "";
    		if(m_useHeaderSAM.getBooleanValue()){
    			params = "-h " + m_headerSAM.getStringValue() + " ";
    		}
     		com = path2samtools + " " + utility + " " + params + path2bamfile + " " + m_inBAM1.getStringValue() + " -o " + outfile;
    	}
    	else if (utility.equals("faidx")){
    		com = path2samtools + " " + utility + " " + path2seqfile;
    	}
    	else if (utility.equals("reheader")) {
    		String outfile = path2bamfile.substring(0,path2bamfile.lastIndexOf(".")) + "_" + utility + ".bam";
       		com = path2samtools + " " + utility + " " + m_rehInSAM.getStringValue() + " " + path2bamfile + " > " + outfile;
    	}
    	else if (utility.equals("merge")) {
    		String outfile = path2bamfile.substring(0,path2bamfile.lastIndexOf("."));
    		outfile += "_" + m_minbam1.getStringValue().substring(m_minbam1.getStringValue().lastIndexOf("/")+1,m_minbam1.getStringValue().lastIndexOf("."));
    		outfile += "_" + utility + ".bam";
    		StringBuffer params = new StringBuffer(6);
    		if(m_mcompression.getBooleanValue()){
    			params.append(" -1 ");    			
    		}
    		if(m_mforce.getBooleanValue()){
    			params.append(" -f ");
    		}
    		if(m_msorted.getBooleanValue()){
    			params.append(" -n ");
    		}
    		if(m_mrgtag.getBooleanValue()){
    			params.append(" -r ");
    		}
    		if(m_muncompressed.getBooleanValue()){
    			params.append(" -u ");
    		}
    		if(m_usemhfile.getBooleanValue()){
    			params.append(" -h "+m_mhfile.getStringValue()+" ");
    		}
    		if(m_usemregion.getBooleanValue()){
    			params.append(" -R "+m_mregion.getStringValue()+" ");
    		}
    		com = path2samtools + " " + utility + " " + params.toString() + " " + outfile + " " + path2bamfile + " " + m_minbam1.getStringValue();
    	}
    	/**else if (utility.equals("faidx")) {
    		com = path2samtools + " " + utility + " " + m_infasta.getStringValue();
    	}**/
    	else if (utility.equals("phase")) {
    		StringBuffer params = new StringBuffer(15);
    		params.append(" -b "+m_prefix.getStringValue()+" ");
    		params.append(" -k "+m_blocklength.getIntValue()+" ");
    		params.append(" -q "+m_hetphred.getIntValue()+" ");
       		params.append(" -Q "+m_minqual.getIntValue()+" ");
       		params.append(" -D "+m_maxdepth.getIntValue()+" ");
       		if(m_fixchimeras.getBooleanValue()){
       			params.append(" -F ");
       		}
       		/*if(m_dropambig.getBooleanValue()){
       			params.append(" -A ");
       		}*/
       		String outfile = path2bamfile.substring(0,path2bamfile.lastIndexOf(".")+1) + "phase";
       		com = "cd "+folder+"; " + path2samtools + " " + utility + " "+ params.toString() + path2bamfile + " >"+outfile; 
    		//redirect STDOUT, else SAM stuff is spammed in log file; it's necessary to "cd" to right directory, else files won't be written there
    	}
   	
    	
    	if(com != "") {
    		// begin QueueSub #################################################
    		if(getAvailableInputFlowVariables().containsKey("JobPrefix")) {
    			String name = getAvailableInputFlowVariables().get("JobPrefix").getStringValue() + "_SamTools";
    			String memory = getAvailableInputFlowVariables().get("Memory").getStringValue();
    			String logfle = path2bamfile.substring(0,path2bamfile.lastIndexOf("/")+1) + "logfile_" + name + ".txt";
    			new QSub(com, name, memory, logfle, true);
     			logBuffer.append("QSub: " + com + "\n");
    			logBuffer.append("See external logfile: " + logfle + "\n");
    		// end QueueSub ###################################################
    		} else {
	    		ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", com);
	        	Process p = b.start();
	        	p.waitFor();
	        	logBuffer.append(ShowOutput.getLogEntry(p, com));
    		}
    	}

    	logBuffer.append(ShowOutput.getNodeEndTime());
    	ShowOutput.writeLogFile(logBuffer);
    	
    	return new BufferedDataTable[]{};
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
    	
    	String utility = m_utility.getStringValue();
    	//flowVariables
        String bamfile = "";
        String fastafile = "";
        String basePathWithFileName = "";
    	
    	//Check OptionalInputPort
		try{
			inSpecs[0].getColumnNames();
			optionalPort=true;
			String[] colnames = inSpecs[0].getColumnNames();

			if(colnames[0].equals("Path2SamTools") && colnames[1].equals("Path2BAMFile")){
			    bamfile = getAvailableInputFlowVariables().get("BAMSAMINFILE").getStringValue();
			    fastafile = getAvailableInputFlowVariables().get("Path2seqFile").getStringValue();
	
				m_samtools.setEnabled(false);
	    		m_bamfile.setEnabled(false);
	    		useSamtools=false;
	    		useBamfile=false;
	    		//useFastafile=false;
	    		if(m_utility.getStringValue().equals("fillmd") || m_utility.getStringValue().equals("faidx")){
	    			m_refseqfile.setEnabled(true);
	    			useFastafile=true;
	    		}


	    		
			}
			else{
				throw new InvalidSettingsException("The in-port is invalid! Don't use an in-port or connect the right node.");
			}

			//m_bfast.setEnabled(false);
    	}catch(NullPointerException npe){
    		
    		
		    bamfile = m_bamfile.getStringValue();
		    fastafile = m_refseqfile.getStringValue();
    		
    		m_samtools.setEnabled(true);
    		m_bamfile.setEnabled(true);
    		useSamtools=true;
    		useBamfile=true;
    		useFastafile=true;
    		if(m_utility.getStringValue().equals("fillmd")){
    				m_refseqfile.setEnabled(true);
    				useFastafile=true;
    		}
			optionalPort=false;
    	
    	}
    	
    	/*******************************************************************
		 ************         Check if Outfile exists            ***********
		 *******************************************************************/
        if(bamfile.length()>0){
        	basePathWithFileName = bamfile.substring(0,bamfile.lastIndexOf("."));
        }
        else if(fastafile.length()>0){
        	basePathWithFileName = fastafile.substring(0,fastafile.lastIndexOf("."));
        }
	    
    	String outfile="";
        if(utility.equals("flagstat") || utility.equals("idxstats") || utility.equals("depth") || utility.equals("rmdup")){
        	outfile = basePathWithFileName + "_" + utility + ".txt";
        }
        else if(utility.equals("faidx")){
        	outfile = fastafile + ".fai";
        	//System.out.println(outfile);
        }
        else if(utility.equals("fillmd")){
        	String outFileExtension = ".sam";
        	if(m_useCompression.getBooleanValue()){
    			outFileExtension = ".bam";
        	}
        	outfile = basePathWithFileName + "_" + utility + outFileExtension;
        }
        else if(utility.equals("fixmate") || utility.equals("rmdup")){
        	outfile = basePathWithFileName + "_" + utility + ".bam";
        }
        else if(utility.equals("merge")){
        	outfile = bamfile.substring(0,bamfile.lastIndexOf("."));
    		outfile += "_" + m_minbam1.getStringValue().substring(m_minbam1.getStringValue().lastIndexOf("/")+1,m_minbam1.getStringValue().lastIndexOf("."));
    		outfile += "_" + utility + ".bam";
        }
        else if(utility.equals("phase")){
        	outfile = bamfile.substring(0,bamfile.lastIndexOf(".")+1) + "phase";
        	String folder = bamfile.substring(0,bamfile.lastIndexOf("/")+1);
        	String outfile_1 = folder + m_prefix.getStringValue() + ".0.bam";
        	String outfile_2 = folder + m_prefix.getStringValue() + ".1.bam";
        	String outfile_3 = folder + m_prefix.getStringValue() + ".chimera.bam";
        	String[] files = {outfile_1,outfile_2,outfile_3};
        	for(int i=0;i<files.length;i++){
        		File o = new File(files[i]);
        		if(o.exists()){
                	setWarningMessage(files[i]+" already exists! Please rename or move to other directory.");
                }
        	}
        }
        else if(utility.equals("reheader")){
        	outfile = bamfile.substring(0,bamfile.lastIndexOf(".")) + "_" + utility + ".bam";
        }


		
        File outpath = new File(outfile);
        if(outpath.exists()){
        	setWarningMessage(outfile+" already exists! Please rename or move to other directory.");
        }
        /***************************************************************************************************/
	
    	
         
    	
     	if(m_refseqfile.getStringValue().length() > 1) {
	     	if(!FileValidator.checkFastaFormat(m_refseqfile.getStringValue())){
	     		throw new InvalidSettingsException("Sequence file is not in fasta format!");
	     	}
     	}
    	

    	/**if(m_utility.getStringValue().equals("faidx") && !FileValidator.checkFastaFormat(m_infasta.getStringValue())){
    		throw new InvalidSettingsException("Error: File not in fasta format!");
    	}**/

        return new DataTableSpec[]{};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {

    	m_utility.saveSettingsTo(settings);
    	m_samtools.saveSettingsTo(settings);
    	m_bamfile.saveSettingsTo(settings);
    	m_refseqfile.saveSettingsTo(settings);
    	//fillmd
    	m_changeIdentBases.saveSettingsTo(settings);
    	m_compression.saveSettingsTo(settings);
    	m_useCompression.saveSettingsTo(settings);
    	m_inputIsSAM.saveSettingsTo(settings);
    	m_modifyQual.saveSettingsTo(settings); 
    	m_bqTag.saveSettingsTo(settings); 
    	m_extendedBAQ.saveSettingsTo(settings); 
    	//m_doCapMapQual.saveSettingsTo(settings); 
    	//m_capMapQual.saveSettingsTo(settings); 
    	//remdup
    	m_removeDup.saveSettingsTo(settings);
    	m_treatPE.saveSettingsTo(settings);
    	//cat
    	m_useHeaderSAM.saveSettingsTo(settings);
    	m_headerSAM.saveSettingsTo(settings);
    	m_inBAM1.saveSettingsTo(settings);
    	//m_inBAM2.saveSettingsTo(settings);
    	//reheader
    	m_rehInSAM.saveSettingsTo(settings);
    	//m_rehInBAM.saveSettingsTo(settings);
    	//merge
    	m_mcompression.saveSettingsTo(settings);
    	m_mforce.saveSettingsTo(settings);
    	m_usemhfile.saveSettingsTo(settings);
    	m_mhfile.saveSettingsTo(settings);
    	m_msorted.saveSettingsTo(settings);
    	m_mregion.saveSettingsTo(settings);
    	m_usemregion.saveSettingsTo(settings);
    	m_mrgtag.saveSettingsTo(settings);
    	m_muncompressed.saveSettingsTo(settings);
    	m_minbam1.saveSettingsTo(settings);
    	//m_minbam2.saveSettingsTo(settings);
    	//faidx
    	//m_infasta.saveSettingsTo(settings);
    	//phase
    	m_blocklength.saveSettingsTo(settings);
    	m_prefix.saveSettingsTo(settings);
    	m_hetphred.saveSettingsTo(settings);
    	m_minqual.saveSettingsTo(settings);
    	m_maxdepth.saveSettingsTo(settings);
    	m_fixchimeras.saveSettingsTo(settings);
    	//m_dropambig.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	
    	m_utility.loadSettingsFrom(settings);
    	m_samtools.loadSettingsFrom(settings);
    	m_bamfile.loadSettingsFrom(settings);
    	m_refseqfile.loadSettingsFrom(settings);
    	//fillmd
    	m_changeIdentBases.loadSettingsFrom(settings);
    	m_compression.loadSettingsFrom(settings);
    	m_useCompression.loadSettingsFrom(settings); 
    	m_inputIsSAM.loadSettingsFrom(settings);
    	m_modifyQual.loadSettingsFrom(settings); 
    	m_bqTag.loadSettingsFrom(settings); 
    	m_extendedBAQ.loadSettingsFrom(settings);
    	//m_doCapMapQual.loadSettingsFrom(settings);
    	//m_capMapQual.loadSettingsFrom(settings);
    	//remdup
    	m_removeDup.loadSettingsFrom(settings);
    	m_treatPE.loadSettingsFrom(settings);
    	//cat
    	m_useHeaderSAM.loadSettingsFrom(settings);
    	m_headerSAM.loadSettingsFrom(settings);
    	m_inBAM1.loadSettingsFrom(settings);
    	//m_inBAM2.loadSettingsFrom(settings);
    	//reheader
    	m_rehInSAM.loadSettingsFrom(settings);
    	//m_rehInBAM.loadSettingsFrom(settings);
    	//merge
    	m_mcompression.loadSettingsFrom(settings);
    	m_mforce.loadSettingsFrom(settings);
    	m_usemhfile.loadSettingsFrom(settings);
    	m_mhfile.loadSettingsFrom(settings);
    	m_msorted.loadSettingsFrom(settings);
    	m_mregion.loadSettingsFrom(settings);
    	m_usemregion.loadSettingsFrom(settings);
    	m_mrgtag.loadSettingsFrom(settings);
    	m_muncompressed.loadSettingsFrom(settings);
    	m_minbam1.loadSettingsFrom(settings);
    	//m_minbam2.loadSettingsFrom(settings);
    	//faidx
    	//m_infasta.loadSettingsFrom(settings);
    	//phase
    	m_blocklength.loadSettingsFrom(settings);
    	m_prefix.loadSettingsFrom(settings);
    	m_hetphred.loadSettingsFrom(settings);
    	m_minqual.loadSettingsFrom(settings);
    	m_maxdepth.loadSettingsFrom(settings);
    	m_fixchimeras.loadSettingsFrom(settings);
    	//m_dropambig.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	
    	m_utility.validateSettings(settings);
    	m_samtools.validateSettings(settings);
    	m_bamfile.validateSettings(settings);
    	m_refseqfile.validateSettings(settings);
    	//fillmd
    	m_changeIdentBases.validateSettings(settings);
    	m_compression.validateSettings(settings);
    	m_useCompression.validateSettings(settings); 
    	m_inputIsSAM.validateSettings(settings);
    	m_modifyQual.validateSettings(settings); 
    	m_bqTag.validateSettings(settings); 
    	m_extendedBAQ.validateSettings(settings); 
      	//m_doCapMapQual.validateSettings(settings); 
    	//m_capMapQual.validateSettings(settings);
    	//remdup
    	m_removeDup.validateSettings(settings);
    	m_treatPE.validateSettings(settings);
    	//cat
    	m_useHeaderSAM.validateSettings(settings);
    	m_headerSAM.validateSettings(settings);
    	m_inBAM1.validateSettings(settings);
    	//m_inBAM2.validateSettings(settings);
    	//reheader
    	m_rehInSAM.validateSettings(settings);
    	//m_rehInBAM.validateSettings(settings);
    	//merge
    	m_mcompression.validateSettings(settings);
    	m_mforce.validateSettings(settings);
    	m_usemhfile.validateSettings(settings);
    	m_mhfile.validateSettings(settings);
    	m_msorted.validateSettings(settings);
    	m_mregion.validateSettings(settings);
    	m_usemregion.validateSettings(settings);
    	m_mrgtag.validateSettings(settings);
    	m_muncompressed.validateSettings(settings);
    	m_minbam1.validateSettings(settings);
    	//m_minbam2.validateSettings(settings);
    	//faidx
    	//m_infasta.validateSettings(settings);
    	//phase
    	m_blocklength.validateSettings(settings);
    	m_prefix.validateSettings(settings);
    	m_hetphred.validateSettings(settings);
    	m_minqual.validateSettings(settings);
    	m_maxdepth.validateSettings(settings);
    	m_fixchimeras.validateSettings(settings);
    	//m_dropambig.validateSettings(settings);
    	
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


package de.helmholtz_muenchen.ibis.ngs.samtools;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;


/**
 * This is the model implementation of SamTools.
 * 
 * @author Sebastian Kopetzky
 * @author Maximilian Hastreiter
 */

public class SamToolsNodeModel extends HTExecutorNodeModel {
    
	
	protected static boolean useSamtools = true;
	protected static boolean useBamfile = true;
	protected static boolean useFastafile = true;
	
	public static final String CFGKEY_UTILITY = "utility";
	public static final String CFGKEY_SAMTOOLS = "samtools";
	public static final String CFGKEY_REFSEQFILE = "refseqfile";
	public static final String CFGKEY_OPTIONAL = "optional";
	
	
	private final SettingsModelString m_utility 			= new SettingsModelString(SamToolsNodeModel.CFGKEY_UTILITY, "");
	private final SettingsModelString m_samtools 			= new SettingsModelString(SamToolsNodeModel.CFGKEY_SAMTOOLS, "");
	private final SettingsModelString m_refseqfile 			= new SettingsModelString(SamToolsNodeModel.CFGKEY_REFSEQFILE, "");
	private final SettingsModelOptionalString m_Optional 	= new SettingsModelOptionalString(CFGKEY_OPTIONAL,"",false);


	/**Parameters for "calmd"**/
//	public static final String CFGKEY_CHANGEIDENTBASES = "changeIdentBases";
//	public static final String CFGKEY_USECOMPRESSION = "useCompression";
//	public static final String CFGKEY_COMPRESSION = "compression";
//	public static final String CFGKEY_INPUTISSAM = "inputIsSam";
//	public static final String CFGKEY_MODIFYQUAL = "modifyQual";
//	public static final String CFGKEY_BQTAG = "bqTag";
//	public static final String CFGKEY_EXTENDEDBAQ = "extendedBAQ";
//	public static final String CFGKEY_DOCAPMAPQUAL = "doCapMapQual";
//	public static final String CFGKEY_CAPMAPQUAL = "capMapQual";
	
//	private final SettingsModelBoolean m_changeIdentBases = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_CHANGEIDENTBASES, false);
//	private final SettingsModelBoolean m_useCompression = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USECOMPRESSION, false);	
//	private final SettingsModelString m_compression = new SettingsModelString(SamToolsNodeModel.CFGKEY_COMPRESSION, "");
//	private final SettingsModelBoolean m_inputIsSAM = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_INPUTISSAM, false);
//	private final SettingsModelBoolean m_modifyQual = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MODIFYQUAL, false);
//	private final SettingsModelBoolean m_bqTag = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_BQTAG, false);
//	private final SettingsModelBoolean m_extendedBAQ = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_EXTENDEDBAQ, false);

	/**Parameters for rmdup**/
	public static final String CFGKEY_REMOVEDUP = "removeDup";
	public static final String CFGKEY_TREATPE = "treatPE";
	private final SettingsModelBoolean m_removeDup = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_REMOVEDUP, false);
	private final SettingsModelBoolean m_treatPE = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_TREATPE, false);
	/**Parameters for cat**/
	public static final String CFGKEY_USEHEADERSAM = "useHeaderSAM";
	public static final String CFGKEY_HEADERSAM = "headerSAM";

	private final SettingsModelBoolean m_useHeaderSAM = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEHEADERSAM, false);
	private final SettingsModelString m_headerSAM = new SettingsModelString(SamToolsNodeModel.CFGKEY_HEADERSAM, "");

	/**Parameters for reheader**/
	public static final String CFGKEY_REHINSAM = "rehInSAM";

	private final SettingsModelString m_rehInSAM = new SettingsModelString(SamToolsNodeModel.CFGKEY_REHINSAM, "");
	
	/**Parameters for merge**/
	public static final String CFGKEY_MCOMPRESSION = "mcompression"; //zlib compression
	public static final String CFGKEY_MFORCE = "mforce";
	public static final String CFGKEY_USEMHFILE = "usemhfile";
	public static final String CFGKEY_MHFILE = "mhfile";
	public static final String CFGKEY_MSORTED = "msorted";
	public static final String CFGKEY_MREGION = "mregion";
	public static final String CFGKEY_USEMREGION = "usemregion";
	public static final String CFGKEY_MRGTAG = "mrgtag";
	public static final String CFGKEY_MUNCOMPRESSED = "muncompressed";
	
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

	/**phase**/
	public static final String CFGKEY_BLOCKLENGTH = "blocklength";
	public static final String CFGKEY_PREFIX = "prefix";
	public static final String CFGKEY_HETPHRED = "hetphred";
	public static final String CFGKEY_MINQUAL = "minqual";
	public static final String CFGKEY_MAXDEPTH = "maxdepth";
	public static final String CFGKEY_FIXCHIMERAS = "fixchimeras";

	private final SettingsModelIntegerBounded m_blocklength = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_BLOCKLENGTH, 13, 1, Integer.MAX_VALUE);
	private final SettingsModelString m_prefix = new SettingsModelString(SamToolsNodeModel.CFGKEY_PREFIX,"phase");
	private final SettingsModelIntegerBounded m_hetphred = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_HETPHRED, 37, 1, Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_minqual = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_MINQUAL, 13, 1, Integer.MAX_VALUE);
	private final SettingsModelIntegerBounded m_maxdepth = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_MAXDEPTH, 256, 1, Integer.MAX_VALUE);
	private final SettingsModelBoolean m_fixchimeras = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_FIXCHIMERAS, false);

	
	protected static Boolean optionalPort 	= false;
	public static String OUT_COL1 			= "Outfile";
	private String OutCellType 				= "FileCell";	
	
	/**
     * Constructor for the node model.
     */
    protected SamToolsNodeModel() {
    
        super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(1));
       
       	addSetting(m_utility);
    	addSetting(m_samtools);
//    	addSetting(m_bamfile);
    	addSetting(m_refseqfile);
    	addSetting(m_Optional);
    	
    	//calmd
//    	addSetting(m_changeIdentBases);
//    	addSetting(m_compression);
//    	addSetting(m_useCompression);
//    	addSetting(m_inputIsSAM);
//    	addSetting(m_modifyQual); 
//    	addSetting(m_bqTag); 
//    	addSetting(m_extendedBAQ); 

    	//remdup
    	addSetting(m_removeDup);
    	addSetting(m_treatPE);
    	//cat
    	addSetting(m_useHeaderSAM);
    	addSetting(m_headerSAM);
//    	addSetting(m_inBAM1);

    	//reheader
    	addSetting(m_rehInSAM);

    	//merge
    	addSetting(m_mcompression);
    	addSetting(m_mforce);
    	addSetting(m_usemhfile);
    	addSetting(m_mhfile);
    	addSetting(m_msorted);
    	addSetting(m_mregion);
    	addSetting(m_usemregion);
    	addSetting(m_mrgtag);
    	addSetting(m_muncompressed);

    	//phase
    	addSetting(m_blocklength);
    	addSetting(m_prefix);
    	addSetting(m_hetphred);
    	addSetting(m_minqual);
    	addSetting(m_maxdepth);
    	addSetting(m_fixchimeras);
    }

    /**
     * {@inheritDoc}
     * @return 
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	    	
    	String path2samtools 	= m_samtools.getStringValue();
    	String path2bamfile 	= "";
    	String path2seqfile 	= m_refseqfile.getStringValue();
    	String lockFile 		= "";
    	String OutCellType 		= "FileCell";		 
    	
    	if(optionalPort){
    		path2bamfile = inData[0].iterator().next().getCell(0).toString();
    	}
    	
    	
    	String basePathWithFileName = "";
    	if(!path2bamfile.isEmpty()){
    	   	basePathWithFileName = path2bamfile.substring(0,path2bamfile.lastIndexOf("."));	
    	}
    	else{
    		basePathWithFileName = path2seqfile.substring(0,path2seqfile.lastIndexOf("."));	
    	}
 
    	String utility = m_utility.getStringValue();
    	ArrayList<String> command = new ArrayList<String>();
    	command.add(path2samtools+" "+utility);
    	
    	if(m_Optional.isActive()){
    		command.add(m_Optional.getStringValue());	
    	}
    	
    	String outfile = "-1";
    	boolean stdOut = false;
    	
    	
    	if(utility.equals("flagstat") || utility.equals("idxstats") || utility.equals("depth")) {
    		command.add(path2bamfile);
    		outfile = basePathWithFileName + "_" + utility + ".txt";
    		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
    		stdOut = true;
    		
//    	} else if(utility.equals("calmd")) {
//
//    		String outFileExtension = ".sam";
//    		OutCellType = "SAMCell";
//    		
//    		if(m_changeIdentBases.getBooleanValue()){
//    			command.add("-e");
//    		}
//    		if(m_useCompression.getBooleanValue()){
//    			outFileExtension = ".bam";
//    			OutCellType = "BAMCell";
//    			if(m_compression.getStringValue().equals("uncompressed")){
//    				command.add("-u");
//    			}
//    			else{
//    				command.add("-b");
//    			}
//    		}
//    		if(m_inputIsSAM.getBooleanValue()){
//    			command.add("-S");
//    		}
//    		if(m_modifyQual.getBooleanValue()){
//    			command.add("-A");
//    		}
//    		if(m_bqTag.getBooleanValue()){
//    			command.add("-r");
//    		}
//    		if(m_extendedBAQ.getBooleanValue()){
//    			command.add("-E");
//    		}
//    		
//    		command.add(path2bamfile);
//    		command.add(path2seqfile);
//    		outfile = basePathWithFileName + "_" + utility + outFileExtension;
//    		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
//    		stdOut = true;
    		
    		
    			
    	} else if (utility.equals("fixmate")) {
    		command.add(path2bamfile);
    		outfile = basePathWithFileName + "_" + utility + ".bam";
    		command.add(outfile);
    		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
    		OutCellType = "BAMCell";
    		
    	} else if (utility.equals("rmdup")) {
    		if(m_removeDup.getBooleanValue()){
    			command.add("-s");
    		}
    		if(m_treatPE.getBooleanValue()){
    			command.add("-S");
    		}
    		command.add(path2bamfile);
    		outfile = basePathWithFileName + "_" + utility + ".bam";
    		command.add(outfile);
    		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
    		OutCellType = "BAMCell";
    		
    	} else if (utility.equals("cat")) {
    		
    		outfile = path2bamfile.substring(0,path2bamfile.lastIndexOf("."));
    		outfile += ".catAll.bam";
    		if(m_useHeaderSAM.getBooleanValue()){
    			command.add("-h " + m_headerSAM.getStringValue());
    		}
    		command.add("-o " + outfile);
    		
            Iterator<DataRow> it = inData[0].iterator();
            int size = 0;
            while(it.hasNext()){
            	command.add(it.next().getCell(0).toString());
            	size++;
            } 
     		if(size<2){
     			throw new InvalidSettingsException("Insufficient number of input files. At least 2 bam files required for samtools cat.");
     		}
     		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
     		OutCellType = "BAMCell";
     		
     		
    	} else if (utility.equals("faidx")){
    		command.add(path2seqfile);
    		outfile = path2seqfile;
    		lockFile = outfile+".faidx" + SuccessfulRunChecker.LOCK_ENDING;
    		
    	} else if (utility.equals("reheader")) {
    		outfile = path2bamfile.substring(0,path2bamfile.lastIndexOf(".")) + "_" + utility + ".bam";
    		command.add(m_rehInSAM.getStringValue());
    		command.add(path2bamfile);
    		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
    		stdOut = true;
    		OutCellType = "BAMCell";
    		
    	} else if (utility.equals("merge")) {
    		outfile = path2bamfile.substring(0,path2bamfile.lastIndexOf("."));
    		outfile += ".mergedAll.bam";
    		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
    		
    		if(m_mcompression.getBooleanValue()){
    			command.add("-1");    			
    		}
    		if(m_mforce.getBooleanValue()){
    			command.add("-f");
    		}
    		if(m_msorted.getBooleanValue()){
    			command.add("-n");
    		}
    		if(m_mrgtag.getBooleanValue()){
    			command.add("-r");
    		}
    		if(m_muncompressed.getBooleanValue()){
    			command.add("-u");
    		}
    		if(m_usemhfile.getBooleanValue()){
    			command.add("-h "+m_mhfile.getStringValue());
    		}
    		if(m_usemregion.getBooleanValue()){
    			command.add("-R "+m_mregion.getStringValue());
    		}
    		
    		command.add(outfile);
    		
            Iterator<DataRow> it = inData[0].iterator();
            int size = 0;
            while(it.hasNext()){
            	command.add(it.next().getCell(0).toString());
            	size++;
            } 
     		if(size<2){
     			throw new InvalidSettingsException("Insufficient number of input files. At least 2 bam files required for samtools cat.");
     		}
     		
    		OutCellType = "BAMCell";
    		
    	} 	else if (utility.equals("phase")) {

    		command.add("-b "+m_prefix.getStringValue());
    		command.add("-k "+m_blocklength.getIntValue());
    		command.add("-q "+m_hetphred.getIntValue());
    		command.add("-Q "+m_minqual.getIntValue());
    		command.add("-D "+m_maxdepth.getIntValue());
       		if(m_fixchimeras.getBooleanValue()){
       			command.add("-F");
       		}

       		outfile = path2bamfile.substring(0,path2bamfile.lastIndexOf(".")+1) + "phase";
       		lockFile = outfile + SuccessfulRunChecker.LOCK_ENDING;
       		command.add(path2bamfile);
       		stdOut = true;
       		OutCellType = "BAMCell";
    	}
   	
    	
    	if(command.size()!=0) {
    		
    		if(!stdOut){	//No extra output file required
    			super.executeCommand(new String[]{StringUtils.join(command, " ")},exec,new File(lockFile),null,outfile+".stdErr");
    			
    		}else{	//StdOut to outfile
    			super.executeCommand(new String[] { StringUtils.join(command, " ") }, exec, null, new File(lockFile), outfile,outfile+".stdErr", null, null, null);
    		}
    	}else{
    		throw new Exception("Something went wrong...Execution command is empty!");
    	}
   	
    	
    	DataColumnSpec dcs = null;
    	if(OutCellType.equals("BAMCell")){
    		dcs = new DataColumnSpecCreator(OUT_COL1, BAMCell.TYPE).createSpec();
    	}else{
    		dcs = new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec();
    	}
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					dcs}));
    	
    	FileCell[] c = new FileCell[]{(FileCell) FileCellFactory.create(outfile)};
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
    	
    	String utility = m_utility.getStringValue();
   	
		if(CompatibilityChecker.inputFileNotOk(m_samtools.getStringValue(), false)) {
			throw new InvalidSettingsException("Set path to samtools binary!");
		}
    	
    	
    	//Check OptionalInputPort
		try{
			inSpecs[0].getColumnNames();
			optionalPort=true;
			
			if(!inSpecs[0].getColumnSpec(0).getType().toString().equals("BAMCell")){
				throw new InvalidSettingsException("Expected BAMCell in first table cell but found a "+inSpecs[0].getColumnSpec(0).getType().toString()+".");
			}
			
			
    	}catch(NullPointerException npe){
    		
    		optionalPort=false;
    		
    		if(!utility.equals("faidx")){
    			setWarningMessage(m_utility.getStringValue()+" requires input table.");
    			throw new InvalidSettingsException(m_utility.getStringValue()+" requires input table.");
    		}	
    	}
    	
    	
     	if(m_refseqfile.getStringValue().length() > 1) {
	     	if(!FileValidator.checkFastaFormat(m_refseqfile.getStringValue())){
	     		throw new InvalidSettingsException("Sequence file is not in fasta format!");
	     	}
     	}

    	DataColumnSpec dcs = null;
    	if(OutCellType.equals("BAMCell")){
    		dcs = new DataColumnSpecCreator(OUT_COL1, BAMCell.TYPE).createSpec();
    	}else{
    		dcs = new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec();
    	}
    	    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					dcs})};
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


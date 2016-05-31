package de.helmholtz_muenchen.ibis.ngs.vep;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of VEP.
 * 
 *
 * @author Tim Jeske
 */
public class VEPNodeModel extends HTExecutorNodeModel {
	
	// the logger instance
    protected static final NodeLogger logger = NodeLogger.getLogger(VEPNodeModel.class);
    
    //annotation
    
	static final String [] TRANSCRIPT_SETS = {"GENCODE Basic","Ensembl","RefSeq"};
	static final String DEF_TRANS_SET = TRANSCRIPT_SETS[0];
	
    static final String CFGKEY_TRANSCRIPT_SET ="transcript_set";
    static final String CFGKEY_CODING_ONLY = "coding_only";
    static final String CFGKEY_FURTHER_OPTIONS = "further_options";
    static final String CFGKEY_SIFT = "sift";
    static final String CFGKEY_POLYPHEN = "polyphen";
    static final String CFGKEY_SYMBOL = "symbol";
    static final String CFGKEY_BIOTYPE = "biotype";
    
	final SettingsModelString m_transcript_set = new SettingsModelString(CFGKEY_TRANSCRIPT_SET,VEPNodeModel.DEF_TRANS_SET);
    final SettingsModelBoolean m_coding_only = new SettingsModelBoolean(CFGKEY_CODING_ONLY,true);
    final SettingsModelBoolean m_sift = new SettingsModelBoolean(VEPNodeModel.CFGKEY_SIFT, true);
    final SettingsModelBoolean m_polyphen = new SettingsModelBoolean(VEPNodeModel.CFGKEY_POLYPHEN, true);
    final SettingsModelBoolean m_symbol = new SettingsModelBoolean(VEPNodeModel.CFGKEY_SYMBOL, true);
    final SettingsModelBoolean m_biotype = new SettingsModelBoolean(VEPNodeModel.CFGKEY_BIOTYPE, true);
    final SettingsModelOptionalString m_further_options = new SettingsModelOptionalString(VEPNodeModel.CFGKEY_FURTHER_OPTIONS, "", false);
    
    //output options
    static final String [] OUT_FORMATS = new String[]{"default","vcf","json","gvf"};
    static final String DEF_FORMAT = OUT_FORMATS[1];
    
    static final String [] STAT_TYPES = {"none", "html", "plain text"};
    static final String DEF_STAT_TYPE = STAT_TYPES[1];
    
    static final String CFGKEY_OUT_FORMAT = "out_format";
    static final String CFGKEY_OUTFOLDER = "outfolder";
    static final String CFGKEY_OVERWRITE = "overwrite";
    static final String CFGKEY_STATS_TYPE = "stats_type";
        
    final SettingsModelString m_out_format = new SettingsModelString(VEPNodeModel.CFGKEY_OUT_FORMAT,VEPNodeModel.DEF_FORMAT);
	final SettingsModelString m_outfolder = new SettingsModelString(CFGKEY_OUTFOLDER,"");
	final SettingsModelBoolean m_overwrite = new SettingsModelBoolean(CFGKEY_OVERWRITE, false);
	final SettingsModelString m_stats_type = new SettingsModelString(CFGKEY_STATS_TYPE,VEPNodeModel.DEF_STAT_TYPE);
	
	//computation options
	static final String DEF_CACHE_DIR = System.getProperty("user.home") + System.getProperty("file.separator") + ".vep";

	static final String CFGKEY_USE_CACHE = "use_cache";
	static final String CFGKEY_CACHE_DIR = "cache_dir";
	static final String CFGKEY_FORKS = "forks";
	static final String CFGKEY_BUFFER = "buffer_size";
	
	final SettingsModelBoolean m_use_cache = new SettingsModelBoolean(CFGKEY_USE_CACHE,true);
	final SettingsModelString m_cache_dir = new SettingsModelString(CFGKEY_CACHE_DIR, DEF_CACHE_DIR);
	final SettingsModelInteger m_forks = new SettingsModelInteger(CFGKEY_FORKS,1);
	final SettingsModelInteger m_buffer = new SettingsModelInteger(VEPNodeModel.CFGKEY_BUFFER, 5000);
	
    //preference page options
    static final String CFGKEY_VEP_PL = "vep_pl";
	final SettingsModelString m_veppl = new SettingsModelString(CFGKEY_VEP_PL,"");
	
	//plugins
	static final String DEF_PLUGIN_DIR = DEF_CACHE_DIR + System.getProperty("file.separator") + "Plugins";
	
	static final String CFGKEY_USE_LOFTEE = "use_loftee";
	final SettingsModelBoolean m_use_loftee = new SettingsModelBoolean(CFGKEY_USE_LOFTEE,false);
	
	static final String CFGKEY_FASTA_FILE = "fasta_file";
	final SettingsModelString m_fasta = new SettingsModelString(CFGKEY_FASTA_FILE,DEF_CACHE_DIR);
	
	static final String CFGKEY_HUMAN_ANCESTOR = "human_ancestor";
	final SettingsModelString m_human_ancestor = new SettingsModelString(CFGKEY_HUMAN_ANCESTOR,"");
	
	static final String CFGKEY_CONSERVATION_FILE = "conservation_file";
	final SettingsModelString m_conservation_file = new SettingsModelString(CFGKEY_CONSERVATION_FILE,"");
	
	static final String CFGKEY_SAMTOOLS_PATH = "samtools_path";
	final SettingsModelString m_samtools_path = new SettingsModelString(CFGKEY_SAMTOOLS_PATH,"");
	
	static final String CFGKEY_FURTHER_PLUGINS = "further_plugins";
	final SettingsModelOptionalString m_further_plugins = new SettingsModelOptionalString(CFGKEY_FURTHER_PLUGINS, "", false);
	
	static final String CFGKEY_PLUGIN_DIR = "plugin_dir";
	final SettingsModelString m_plugin_dir = new SettingsModelString(CFGKEY_PLUGIN_DIR, DEF_PLUGIN_DIR);
	
	public static final String OUT_COL1 = "Path2VEP_AnnotationVCF";
	
	private int vcf_index;
	private DataType outType = FileCell.TYPE;
	
    /**
     * Constructor for the node model.
     */
    protected VEPNodeModel() {
    	super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
    	
    	boolean use_hte = IBISKNIMENodesPlugin.getBooleanPreference(IBISKNIMENodesPlugin.USE_HTE);
    	
    	if(use_hte) {
    		m_overwrite.setBooleanValue(true);
    		m_overwrite.setEnabled(false);
    	}
    	
    	addSetting(m_out_format);
    	addSetting(m_veppl);
        addSetting(m_fasta);
        addSetting(m_outfolder);
        addSetting(m_stats_type);
        addSetting(m_coding_only);
        addSetting(m_further_options);
        addSetting(m_use_cache);
        addSetting(m_cache_dir);
        addSetting(m_plugin_dir);
        addSetting(m_further_plugins);
        addSetting(m_use_loftee);
        addSetting(m_human_ancestor);
        addSetting(m_conservation_file);
        addSetting(m_samtools_path);
        addSetting(m_forks);
        addSetting(m_transcript_set);
        addSetting(m_overwrite);
        addSetting(m_buffer);
        addSetting(m_sift);
        addSetting(m_polyphen);
        addSetting(m_symbol);
        addSetting(m_biotype);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	//check input file
    	String vcf_infile = inData[0].iterator().next().getCell(vcf_index).toString();
    	
    	if(CompatibilityChecker.inputFileNotOk(vcf_infile)) {
    		throw new InvalidSettingsException("Input VCF file does not exist!");
    	}
    	
    	String outfileBase = "";
    	String outfolder = m_outfolder.getStringValue();
    	if(outfolder.equals("") || Files.notExists(Paths.get(outfolder))) {
    		outfileBase = vcf_infile;
    	} else {
    		outfileBase = m_outfolder.getStringValue()+ System.getProperty("file.separator")+ new File(vcf_infile).getName();
    	}
    	outfileBase = outfileBase.split("\\.vcf")[0];
    	
    	
    	/**create command**/
    	String cmd = "perl "+m_veppl.getStringValue();
    	cmd += " -i " + vcf_infile;
    	
    	//output file
    	String out_type = m_out_format.getStringValue();
    	String res_file = outfileBase + ".vep." + out_type;
    	cmd += " -o " + res_file;
    	if(!out_type.equals("default")) {
    		cmd += " --" + out_type;
    	}
    	
    	//stats file
    	String stats_type = m_stats_type.getStringValue();
    	String stats_file;
    	if(stats_type.equals("html")) {
    		stats_file = outfileBase + ".vep_stats.html";
    		cmd += " --stats_file " + stats_file;
    	} else if(stats_type.equals("plain text")) {
    		stats_file = outfileBase + ".vep_stats.txt";
    		cmd += " --stats_text " + stats_file;
    	} else {
    		cmd += " --no_stats";
    	}
    	
    	String transcriptSet= m_transcript_set.getStringValue();
    	if(transcriptSet.equals("RefSeq")) {
    		cmd += " --refseq";
    	} else if(transcriptSet.equals("GENCODE Basic")) {
    		cmd += " --gencode_basic";
    	}
    	
    	if(m_coding_only.getBooleanValue()) {
    		cmd += " --coding_only";
    	}
    	
    	if(m_forks.getIntValue()>1) {
    		cmd += " --fork " + m_forks.getIntValue();
    	}
    	cmd += " --buffer_size " + m_buffer.getIntValue();
    	
    	if(m_overwrite.getBooleanValue()) {
    		cmd += " --force_overwrite";
    	}
    	
    	//cache usage
    	String cache_dir = m_cache_dir.getStringValue();
    	if(m_use_cache.getBooleanValue()) {
    		if(cache_dir.equals("") || Files.notExists(Paths.get(cache_dir))) {
    			setWarningMessage("Cache directory was not specified!");
    		} else {
    			cmd += " --dir_cache " + cache_dir; 
    		}
    		cmd += " --offline";
    	} else {
    		cmd += " --database";
    	}
    	
    	
    	//default VEP parameters
    	cmd += " --no_progress";
    	cmd += " --allele_number";
    	
    	if(m_sift.getBooleanValue()) {
    		cmd += " --sift b";
    	}
    	
    	if(m_polyphen.getBooleanValue()) {
    		cmd += " --polyphen b";
    	}
    	
    	if(m_symbol.getBooleanValue()) {
        	cmd += " --symbol";
    	}
    	
    	if(m_biotype.getBooleanValue()) {
        	cmd += " --biotype";
    	}
    	
    	if(m_further_options.isActive()) {
    		cmd += " " + m_further_options.getStringValue().trim();
    	}    	
    	
    	//fasta file
    	String fasta_file = m_fasta.getStringValue();
    	String [] flagsRequiringFastAFile = {"--hgvs", "--check_ref"};
    	boolean requireFastA = false;
    	for(String s: flagsRequiringFastAFile) {
    		if(cmd.contains(s)) {
    			requireFastA = true;
    			break;
    		}
    	}
    	if(CompatibilityChecker.inputFileNotOk(fasta_file) && (requireFastA  || m_use_loftee.getBooleanValue())) {
			throw new InvalidSettingsException("Chosen settings require FastA file for looking up reference sequence!");
		} 
    	
    	if(CompatibilityChecker.inputFileNotOk(fasta_file)) {
    		setWarningMessage("Given path to FastA file invalid!");
    	} else {
			cmd += " --fasta "+fasta_file;
		}
    	
    	//plugin parameters
    	String plugin_dir = m_plugin_dir.getStringValue();
    	if(plugin_dir.equals("")|| Files.notExists(Paths.get(plugin_dir))) {
    		setWarningMessage("Plugin directory was not specified!");
    	} else {
    		cmd += " --dir_plugins " + plugin_dir;
    	}
    	
    	String path_variable = "PATH="+System.getenv("PATH") ;
    	String samtools_path = m_samtools_path.getStringValue();
    	if(samtools_path.equals("")|| Files.notExists(Paths.get(samtools_path))) {
    		setWarningMessage("Samtools PATH was not specified!");
    	} else {
    		//remove samtools
    		int pos = samtools_path.lastIndexOf(System.getProperty("file.separator"));
    		path_variable +=":"+samtools_path.substring(0,pos);
    	}
    	
    	//LOFTEE parameters

    	if(m_use_loftee.getBooleanValue()) {
    		cmd += " --plugin LoF";
    		String human_ancestor = m_human_ancestor.getStringValue();
    		if(human_ancestor.equals("")|| Files.notExists(Paths.get(human_ancestor))) {
    			setWarningMessage("Human ancestor file was not specified!");
    		} else {
    			cmd += ",human_ancestor_fa:"+human_ancestor;
    		}
    		if(m_forks.getIntValue()<=1) {
    			String conservation_file = m_conservation_file.getStringValue();
    			if(conservation_file.equals("")|| Files.notExists(Paths.get(conservation_file))) {
    				setWarningMessage("Phylocsf.sql was not specified!");
    			} else {
    				cmd += ",conservation_file:"+conservation_file;
    			}
    		}
    	}
    	
    	//user defined plugins
    	if(m_further_plugins.isActive()) {
    		cmd += " " + m_further_plugins.getStringValue().trim();
    	}
    	
    	String stdOutFile = outfileBase + ".vep.stdout";
    	
    	String stdErrFile = outfileBase + ".vep.stderr";
    	
    	File lockFile = new File(outfileBase + ".vep" + SuccessfulRunChecker.LOCK_ENDING);
    	
    	String perl5lib_variable = "PERL5LIB="+System.getenv("PERL5LIB");

    	logger.info("COMMAND: "+cmd);
    	logger.info("ENVIRONMENT: "+path_variable);
    	logger.info("PERL5LIB: "+perl5lib_variable);
    	super.executeCommand(new String[]{cmd}, exec, new String[]{path_variable, perl5lib_variable}, lockFile, stdOutFile, stdErrFile, null, null, null);
    	
    	
    	//Create Output Table
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, outType).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(res_file)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	
        return new BufferedDataTable[]{outTable};
    }



    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	String script_warning = CheckUtils.checkSourceFile(m_veppl.getStringValue());
    	if(script_warning != null) {
    		setWarningMessage(script_warning);
    	}
    	
//    	String outfolder_warning = CheckUtils.checkDestinationDirectory(m_outfolder.getStringValue());
//    	if(outfolder_warning!=null) {
//    		setWarningMessage(outfolder_warning);
//    	}
    	
    	if(m_use_loftee.getBooleanValue()) {
    		try {
	    		 //Version control
	            if(FileValidator.versionControl(m_samtools_path.getStringValue(),"SAMTOOLS")==1){
	            	setWarningMessage("WARNING: You are using a newer SAMTOOLS version than "+FileValidator.SAMTOOLS_VERSION +"! This may cause problems!");
	            }else if(FileValidator.versionControl(m_samtools_path.getStringValue(),"SAMTOOLS")==2){
	            	setWarningMessage("WARNING: You are using an older SAMTOOLS version than "+FileValidator.SAMTOOLS_VERSION +"! This may cause problems!");
	            }else if(FileValidator.versionControl(m_samtools_path.getStringValue(),"SAMTOOLS")==-1){
	            	setWarningMessage("Your samtools version could not be determined! Correct behaviour can only be ensured for samtools version "+FileValidator.SAMTOOLS_VERSION+".");
	            }
    		} catch (Exception e) {
    			throw new InvalidSettingsException("Specify a valid SAMTOOLS version!");
    		}
    	}
    	
    	vcf_index = -1;
    	for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
    		if(inSpecs[0].getColumnSpec(i).getType().toString().equals("VCFCell")) {
    			vcf_index = i;
    		}
    	}
    	
    	if(vcf_index==-1) {
    		throw new InvalidSettingsException("This node is not compatible with the precedent node as there is no VCF file in the input table!");
    	}
    	
    	if(m_out_format.getStringValue().equals("vcf")) {
    		outType = VCFCell.TYPE;
    	}
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, outType).createSpec()})};
    }
}
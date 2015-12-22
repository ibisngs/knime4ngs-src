package de.helmholtz_muenchen.ibis.ngs.vep;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

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
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.CheckUtils;


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
    
    //general options
	static final String DEF_CACHE_DIR = System.getProperty("user.home") + System.getProperty("file.separator") + ".vep";
	static final String DEF_PLUGIN_DIR = DEF_CACHE_DIR + System.getProperty("file.separator") + "Plugins";

    
    static final String CFGKEY_VEP_PL = "vep_pl";
	final SettingsModelString m_veppl = new SettingsModelString(CFGKEY_VEP_PL,"");
	
	static final String CFGKEY_FASTA_FILE = "fasta_file";
	final SettingsModelString m_fasta = new SettingsModelString(CFGKEY_FASTA_FILE,DEF_CACHE_DIR);
	
	static final String CFGKEY_OUTFOLDER = "outfolder";
	final SettingsModelString m_outfolder = new SettingsModelString(CFGKEY_OUTFOLDER,"");
	
	static final String CFGKEY_OVERWRITE = "overwrite";
	private final SettingsModelBoolean m_overwrite = new SettingsModelBoolean(CFGKEY_OVERWRITE, false);
	
	static final String CFGKEY_USE_CACHE = "use_cache";
	final SettingsModelBoolean m_use_cache = new SettingsModelBoolean(CFGKEY_USE_CACHE,true);
	
	static final String CFGKEY_FORKS = "forks";
	private final SettingsModelInteger m_forks = new SettingsModelInteger(CFGKEY_FORKS,1);
	
	static final String CFGKEY_CODING_ONLY = "coding_only";
	final SettingsModelBoolean m_coding_only = new SettingsModelBoolean(CFGKEY_CODING_ONLY,true);
	
	static final String CFGKEY_TRANSCRIPT_SET ="transcript_set";
	static final String [] TRANSCRIPT_SETS = {"GENCODE Basic","Ensembl","RefSeq"};
	final SettingsModelString m_transcript_set = new SettingsModelString(CFGKEY_TRANSCRIPT_SET,"");
	
	//advanced options: cache directory, plugin directory, stats type
	static final String CFGKEY_STATS_TYPE = "stats_type";
	static final String [] STAT_TYPES = {"none", "html", "plain text"};
	final SettingsModelString m_stats_type = new SettingsModelString(CFGKEY_STATS_TYPE,"");
	
	static final String CFGKEY_CACHE_DIR = "cache_dir";
	final SettingsModelString m_cache_dir = new SettingsModelString(CFGKEY_CACHE_DIR, DEF_CACHE_DIR);
	
	static final String CFGKEY_PLUGIN_DIR = "plugin_dir";
	final SettingsModelString m_plugin_dir = new SettingsModelString(CFGKEY_PLUGIN_DIR, DEF_PLUGIN_DIR);
	
	//loftee
	static final String CFGKEY_USE_LOFTEE = "use_loftee";
	final SettingsModelBoolean m_use_loftee = new SettingsModelBoolean(CFGKEY_USE_LOFTEE,false);
	
	static final String CFGKEY_HUMAN_ANCESTOR = "human_ancestor";
	final SettingsModelString m_human_ancestor = new SettingsModelString(CFGKEY_HUMAN_ANCESTOR,"");
	
	static final String CFGKEY_CONSERVATION_FILE = "conservation_file";
	final SettingsModelString m_conservation_file = new SettingsModelString(CFGKEY_CONSERVATION_FILE,"");
	
	static final String CFGKEY_SAMTOOLS_PATH = "samtools_path";
	final SettingsModelString m_samtools_path = new SettingsModelString(CFGKEY_SAMTOOLS_PATH,"");
	
	public static final String OUT_COL1 = "Path2VEP_AnnotationVCF";
	
	private int vcf_index;
	
    /**
     * Constructor for the node model.
     */
    protected VEPNodeModel() {
    	super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	//perl script
    	String script = m_veppl.getStringValue();
    	
    	//input file
    	String vcf_infile = inData[0].iterator().next().getCell(vcf_index).toString();
    	
    	if(Files.notExists(Paths.get(vcf_infile))) {
    		throw new InvalidSettingsException("Input VCF file does not exist!");
    	}
    	
    	String outfileBase = m_outfolder.getStringValue()+ System.getProperty("file.separator")+ new File(vcf_infile).getName();
    	outfileBase = outfileBase.split("\\.vcf")[0];
    	
    	//output file + stats file
    	String res_file,stats_file,stats_type;
    	res_file = outfileBase + ".vep.vcf";
    	
    	
    	String cmd = "perl "+script;
    	cmd += " -i " + vcf_infile;
    	cmd += " -o " + res_file;
    		
    	stats_type = m_stats_type.getStringValue();
    	if(stats_type.equals("html")) {
    		stats_file = outfileBase + ".vep_stats.html";
    		cmd += " --stats_file " + stats_file;
    	} else if(stats_type.equals("plain text")) {
    		stats_file = outfileBase + ".vep_stats.txt";
    		cmd += " --stats_text " + stats_file;
    	} else {
    		stats_file = "";
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
    	
    	//fasta file
    	String fasta_file = m_fasta.getStringValue();
    	if(fasta_file.equals("")|| Files.notExists(Paths.get(fasta_file))) {
			setWarningMessage("No fasta file specified for looking up reference sequence!");
		} else {
			cmd += " --fasta "+fasta_file;
		}
    	
    	//default VEP parameters
    	cmd += " --no_progress";
    	cmd += " --vcf";
    	cmd += " --buffer_size 10000";
    	cmd += " --allele_number";
    	cmd += " --sift b";
    	cmd += " --polyphen b";
    	cmd += " --symbol";
    	cmd += " --numbers";
    	cmd += " --biotype";
    	cmd += " --total_length";
    	cmd += " --canonical";
    	cmd += " --ccds";
    	cmd += " --fields Consequence,ALLELE_NUM,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,CCDS,LoF,LoF_filter,LoF_flags";
    	
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
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()}));
    	
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
    	
    	String outfolder_warning = CheckUtils.checkDestinationDirectory(m_outfolder.getStringValue());
    	if(outfolder_warning!=null) {
    		setWarningMessage(outfolder_warning);
    	}
    	
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
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	super.saveSettingsTo(settings);
        m_veppl.saveSettingsTo(settings);
        m_fasta.saveSettingsTo(settings);
        m_outfolder.saveSettingsTo(settings);
        m_stats_type.saveSettingsTo(settings);
        m_coding_only.saveSettingsTo(settings);
        m_use_cache.saveSettingsTo(settings);
        m_cache_dir.saveSettingsTo(settings);
        m_plugin_dir.saveSettingsTo(settings);
        m_use_loftee.saveSettingsTo(settings);
        m_human_ancestor.saveSettingsTo(settings);
        m_conservation_file.saveSettingsTo(settings);
        m_samtools_path.saveSettingsTo(settings);
        m_forks.saveSettingsTo(settings);
        m_transcript_set.saveSettingsTo(settings);
        m_overwrite.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
        m_veppl.loadSettingsFrom(settings);
        m_fasta.loadSettingsFrom(settings);
        m_outfolder.loadSettingsFrom(settings);
        m_stats_type.loadSettingsFrom(settings);
        m_coding_only.loadSettingsFrom(settings);
        m_use_cache.loadSettingsFrom(settings);
        m_cache_dir.loadSettingsFrom(settings);
        m_plugin_dir.loadSettingsFrom(settings);
        m_use_loftee.loadSettingsFrom(settings);
        m_human_ancestor.loadSettingsFrom(settings);
        m_conservation_file.loadSettingsFrom(settings);
        m_samtools_path.loadSettingsFrom(settings);
        m_forks.loadSettingsFrom(settings);
        m_transcript_set.loadSettingsFrom(settings);
        m_overwrite.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.validateSettings(settings);
        m_veppl.validateSettings(settings);
        m_fasta.validateSettings(settings);
        m_outfolder.validateSettings(settings);
        m_stats_type.validateSettings(settings);
        m_coding_only.validateSettings(settings);
        m_use_cache.validateSettings(settings);
        m_cache_dir.validateSettings(settings);
        m_plugin_dir.validateSettings(settings);
        m_use_loftee.validateSettings(settings);
        m_human_ancestor.validateSettings(settings);
        m_conservation_file.validateSettings(settings);
        m_samtools_path.validateSettings(settings);
        m_forks.validateSettings(settings);
        m_transcript_set.validateSettings(settings);
        m_overwrite.validateSettings(settings);
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

	@Override
	protected void reset() {
		
	}

}


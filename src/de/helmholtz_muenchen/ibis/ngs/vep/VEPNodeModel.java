package de.helmholtz_muenchen.ibis.ngs.vep;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

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
import org.knime.core.node.defaultnodesettings.SettingsModel;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.ngs.vat.VATNodeModel;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of VEP.
 * 
 *
 * @author tim.jeske
 */
public class VEPNodeModel extends HTExecutorNodeModel {
	
	// the logger instance
    protected static final NodeLogger logger = NodeLogger.getLogger(VATNodeModel.class);
    
    ArrayList<SettingsModel> models = new ArrayList<SettingsModel>();
    
    //general options
    static final String CFGKEY_VEP_PL = "vep_pl";
	final SettingsModelString m_veppl = new SettingsModelString(CFGKEY_VEP_PL,"");
	
    static final String CFGKEY_VCF_INFILE = "vcf_infile";
	final SettingsModelString m_vcfin = new SettingsModelString(CFGKEY_VCF_INFILE,"");
	
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
	static final String [] TRANSCRIPT_SETS = {"Ensembl","RefSeq","GENCODE Basic"};
	final SettingsModelString m_transcript_set = new SettingsModelString(CFGKEY_TRANSCRIPT_SET,"");
	
	//advanced options: cache directory, plugin directory, stats type
	static final String CFGKEY_STATS_TYPE = "stats_type";
	static final String [] STAT_TYPES = {"none", "html", "plain text"};
	final SettingsModelString m_stats_type = new SettingsModelString(CFGKEY_STATS_TYPE,"");
	
	static final String CFGKEY_CACHE_DIR = "cache_dir";
	static final String DEF_CACHE_DIR = System.getProperty("user.home") + System.getProperty("file.separator") + ".vep";
	final SettingsModelString m_cache_dir = new SettingsModelString(CFGKEY_CACHE_DIR, DEF_CACHE_DIR);
	
	static final String CFGKEY_PLUGIN_DIR = "plugin_dir";
	static final String DEF_PLUGIN_DIR = DEF_CACHE_DIR + System.getProperty("file.separator") + "Plugins";
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
	
	
	//tabix path
	static final String CFGKEY_TABIX_PATH = "tabix_path";
	final SettingsModelString m_tabix_path = new SettingsModelString(CFGKEY_TABIX_PATH,"");
	
	//CADD
	static final String CFGKEY_USE_CADD = "use_cadd";
	final SettingsModelBoolean m_use_cadd = new SettingsModelBoolean(CFGKEY_USE_CADD,false);
	
	static final String CFGKEY_FIRST_CADD_FILE = "first_cadd_file";
	final SettingsModelString m_first_cadd_file = new SettingsModelString(CFGKEY_FIRST_CADD_FILE,"");
	
	static final String CFGKEY_SEC_CADD_FILE = "sec_cadd_file";
	final SettingsModelString m_sec_cadd_file = new SettingsModelString(CFGKEY_SEC_CADD_FILE,"");
	
	//ExAC
	static final String CFGKEY_USE_EXAC = "use_exac";
	final SettingsModelBoolean m_use_exac = new SettingsModelBoolean(CFGKEY_USE_EXAC,false);
	
	static final String CFGKEY_EXAC_FILE = "exac_file";
	final SettingsModelString m_exac_file = new SettingsModelString(CFGKEY_EXAC_FILE,"");
	
	public static final String OUT_COL1 = "Path2VEP_AnnotationVCF";
	public static final String OUT_COL2 = "Path2VEP_StatsFile";
	
	public static boolean optionalPort=false;
	
    /**
     * Constructor for the node model.
     */
    protected VEPNodeModel() {
    	super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(1));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	String cmd;
    	
    	
    	//perl script
    	String script = m_veppl.getStringValue();
    	if(script.equals("") || Files.notExists(Paths.get(script))) {
    		System.err.println("Path to variant_effect_predictor.pl was not specified!");
    	}
    	cmd = "perl "+script;
    	
    	//input file
    	String vcf_infile;
    	if(optionalPort){	//Input Table available
    		//Get File via Table
    		vcf_infile = inData[0].iterator().next().getCell(0).toString();
    	}else{
    		//Get File via FileSelector
    		vcf_infile = m_vcfin.getStringValue();
    		if(vcf_infile.equals("") || Files.notExists(Paths.get(vcf_infile))) {
    			System.err.println("No input vcf file specified!");
    		}
    	}
    	cmd += " -i " + vcf_infile;
    	
    	String outfileBase = m_outfolder.getStringValue()+ System.getProperty("file.separator")+ new File(vcf_infile).getName();
    	
    	//output file + stats file
    	String res_file,stats_file,stats_type;
    	res_file = outfileBase.replace("vcf", "vep.vcf");
    	cmd += " -o " + res_file;
    		
    	stats_type = m_stats_type.getStringValue();
    	if(stats_type.equals("html")) {
    		stats_file = outfileBase.replace("vcf", "vep_stats.html");
    		cmd += " --stats_file " + stats_file;
    	} else if(stats_type.equals("plain text")) {
    		stats_file = outfileBase.replace("vcf", "vep_stats.txt");
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
    			System.err.println("Cache directory was not specified!");
    		} else {
    			cmd += " --dir_cache " + cache_dir; 
    		}
    		cmd += " --offline";
    	} else {
    		cmd += " --database";
    	}
    	
    	//default VEP parameters
    	cmd += " --no_progress";
    	cmd += " --vcf";
    	cmd += " --pubmed";
    	
    	
    	//plugin parameters
    	String plugin_dir = m_plugin_dir.getStringValue();
    	if(plugin_dir.equals("")|| Files.notExists(Paths.get(plugin_dir))) {
    		System.err.println("Plugin directory was not specified!");
    	} else {
    		cmd += " --dir_plugins " + plugin_dir;
    	}
    	
    	String environment = "PATH="+System.getenv("PATH") ;
    	String samtools_path = m_samtools_path.getStringValue();
    	if(samtools_path.equals("")|| Files.notExists(Paths.get(samtools_path))) {
    		System.err.println("Samtools PATH was not specified!");
    	} else {
    		environment +=":"+samtools_path;
    	}
    	String tabix_path = m_tabix_path.getStringValue();
    	if(tabix_path.equals("")|| Files.notExists(Paths.get(tabix_path))) {
    		System.err.println("Tabix path was not specified!");
    	} else {
    		environment += ":"+tabix_path;
    	}
    	//LOFTEE parameters

    	if(m_use_loftee.getBooleanValue()) {
    		cmd += " --plugin LoF";
    		String human_ancestor = m_human_ancestor.getStringValue();
    		if(human_ancestor.equals("")|| Files.notExists(Paths.get(human_ancestor))) {
    			System.err.println("Human ancestor file was not specified!");
    		} else {
    			cmd += ",human_ancestor_fa:"+human_ancestor;
    		}
    		if(m_forks.getIntValue()<=1) {
    			String conservation_file = m_conservation_file.getStringValue();
    			if(conservation_file.equals("")|| Files.notExists(Paths.get(conservation_file))) {
    				System.err.println("Phylocsf.sql was not specified!");
    			} else {
    				cmd += ",conservation_file:"+conservation_file;
    			}
    		}
    	}
    	
    	//CADD parameters
    	String first_cadd_file = m_first_cadd_file.getStringValue();
    	String sec_cadd_file = m_sec_cadd_file.getStringValue();
    	if(first_cadd_file.equals("") && sec_cadd_file.equals("")) {
    		System.err.println("No CADD file specified!");
    	}
    	if(m_use_cadd.getBooleanValue()) {
    		cmd += " --plugin CADD";
    		if(!first_cadd_file.equals("") && Files.exists(Paths.get(first_cadd_file))) {
    			cmd += ","+first_cadd_file;
    		}
    		if(!sec_cadd_file.equals("") && Files.exists(Paths.get(sec_cadd_file))) {
    			cmd += ","+sec_cadd_file;
    		}
    	}
    	
    	//ExAC parameters
    	String exac_file = m_exac_file.getStringValue();
    	if(m_use_exac.getBooleanValue()) {
    		cmd += " --plugin ExAC";
    	}
    	if(exac_file.equals("")|| Files.notExists(Paths.get(exac_file))) {
			System.err.println("No ExAC file specified!");
		} else {
			cmd += ","+exac_file;
		}
    	
    	
    	String stdOutFile = outfileBase.replace("vcf", "vep.stdout");
    	String stdErrFile = outfileBase.replace("vcf", "vep.stderr");
    	File lockFile = new File(outfileBase.replace("vcf", "vep" + SuccessfulRunChecker.LOCK_ENDING));

    	System.out.println("COMMAND: "+cmd);
    	System.out.println("ENVIRONMENT: "+environment);
    	
    	super.executeCommand(new String[]{cmd}, exec, new String[]{environment}, logger, lockFile, stdOutFile, stdErrFile, null, null, null);
    	
    	
    	//Create Output Table
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(res_file),
    			(FileCell) FileCellFactory.create(stats_file)};
    	
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
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

    	try{
			inSpecs[0].getColumnNames();
			optionalPort=true;
			
		}catch(NullPointerException e){}

        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	super.saveSettingsTo(settings);
        m_veppl.saveSettingsTo(settings);
        m_vcfin.saveSettingsTo(settings);
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
        m_tabix_path.saveSettingsTo(settings);
        m_use_cadd.saveSettingsTo(settings);
        m_first_cadd_file.saveSettingsTo(settings);
        m_sec_cadd_file.saveSettingsTo(settings);
        m_use_exac.saveSettingsTo(settings);
        m_exac_file.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
        m_veppl.loadSettingsFrom(settings);
        m_vcfin.loadSettingsFrom(settings);
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
        m_tabix_path.loadSettingsFrom(settings);
        m_use_cadd.loadSettingsFrom(settings);
        m_first_cadd_file.loadSettingsFrom(settings);
        m_sec_cadd_file.loadSettingsFrom(settings);
        m_use_exac.loadSettingsFrom(settings);
        m_exac_file.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.validateSettings(settings);
        m_veppl.validateSettings(settings);
        m_vcfin.validateSettings(settings);
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
        m_tabix_path.validateSettings(settings);
        m_use_cadd.validateSettings(settings);
        m_first_cadd_file.validateSettings(settings);
        m_sec_cadd_file.validateSettings(settings);
        m_use_exac.validateSettings(settings);
        m_exac_file.validateSettings(settings);
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


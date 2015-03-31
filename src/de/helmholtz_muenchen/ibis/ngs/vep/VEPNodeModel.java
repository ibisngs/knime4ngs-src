package de.helmholtz_muenchen.ibis.ngs.vep;

import java.io.File;
import java.io.IOException;

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
    
	static final String CFGKEY_VCF_INFILE = "vcf_infile";
	final SettingsModelString m_vcfin = new SettingsModelString(CFGKEY_VCF_INFILE,"");
	
	static final String CFGKEY_VEP_PL = "vep_pl";
	final SettingsModelString m_veppl = new SettingsModelString(CFGKEY_VEP_PL,"");
	
	static final String CFGKEY_STATS_TYPE = "stats_type";
	static final String [] STAT_TYPES = {"none", "html", "plain text"};
	final SettingsModelString m_stats_type = new SettingsModelString(CFGKEY_STATS_TYPE,"");
	
	static final String CFGKEY_CODING_ONLY = "coding_only";
	final SettingsModelBoolean m_coding_only = new SettingsModelBoolean(CFGKEY_CODING_ONLY,true);
	
	static final String CFGKEY_FORKS = "forks";
	private final SettingsModelInteger m_forks = new SettingsModelInteger(CFGKEY_FORKS,1);
	
	//cache usage
	static final String CFGKEY_USE_CACHE = "use_cache";
	final SettingsModelBoolean m_use_cache = new SettingsModelBoolean(CFGKEY_USE_CACHE,true);
	
	static final String CFGKEY_CACHE_DIR = "cache_dir";
	static final String DEF_CACHE_DIR = System.getProperty("user.home") + System.getProperty("file.separator") + ".vep";
	final SettingsModelString m_cache_dir = new SettingsModelString(CFGKEY_CACHE_DIR, DEF_CACHE_DIR);
	
	//pugin directory
	static final String CFGKEY_PLUGIN_DIR = "plugin_dir";
	final SettingsModelString m_plugin_dir = new SettingsModelString(CFGKEY_PLUGIN_DIR, DEF_CACHE_DIR);
	
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

    	String vcf_infile, res_file, stats_file, stdOutFile, stdErrFile, cmd;
    	String environment = "";
    	File lockFile;
    	
    	//perl script
    	cmd = "perl "+m_veppl.getStringValue();
    	
    	//input file
    	if(optionalPort){	//Input Table available
    		//Get File via Table
    		vcf_infile = inData[0].iterator().next().getCell(0).toString();
    	}else{
    		//Get File via FileSelector
    		vcf_infile = m_vcfin.getStringValue();
    	}
    	cmd += " -i " + vcf_infile;
    	
    	//output file + stats file
    	res_file = vcf_infile.replace("vcf", "VEP_annotation.vcf");
    	cmd += " -o " + res_file;
    		
    	if(m_stats_type.getStringValue().equals("html")) {
    		stats_file = vcf_infile.replace("vcf", "VEP_stats.html");
    		cmd += " --stats_file " + stats_file;
    	} else if(m_stats_type.getStringValue().equals("plain text")) {
    		stats_file = vcf_infile.replace("vcf", "VEP_stats.txt");
    		cmd += " --stats_text " + stats_file;
    	} else {
    		stats_file = "";
    		cmd += " --no_stats";
    	}
    	
    	if(m_coding_only.getBooleanValue()) {
    		cmd += " --coding_only";
    	}
    	
    	if(m_forks.getIntValue()>1) {
    		cmd += " --fork " + m_forks.getIntValue();
    	}
    	
    	//cache usage
    	if(m_use_cache.getBooleanValue()) {
    		cmd += " --dir_cache " + m_cache_dir.getStringValue(); 
    		cmd += " --offline";
    	} else {
    		cmd += " --database";
    	}
    	
    	//default VEP parameters
    	cmd += " --no_progress";
    	cmd += " --vcf";
    	
    	//LOFTEE parameters
    	if(m_use_loftee.getBooleanValue()) {
    		cmd += " --dir_plugins " + m_plugin_dir.getStringValue();
    		cmd += " --plugin LoF";
    		cmd += ",human_ancestor_fa:"+m_human_ancestor.getStringValue();
    		if(m_forks.getIntValue()<=1) {
    			cmd += ",conservation_file:"+m_conservation_file.getStringValue();
    		}
    		
    		environment = "PATH=$PATH:"+m_samtools_path.getStringValue();
    	}
    	
    	stdOutFile = vcf_infile.replace("vcf","VEP_stdout.vcf");
    	stdErrFile = vcf_infile.replace("vcf", "VEP_stderr.vcf");
    	lockFile = new File(vcf_infile.substring(0,vcf_infile.lastIndexOf(".")) + ".VEP" +  SuccessfulRunChecker.LOCK_ENDING);

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


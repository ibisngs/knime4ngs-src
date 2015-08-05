package de.helmholtz_muenchen.ibis.ngs.vcffilter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;

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
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;
/**
 * This is the model implementation of LOFFilter.
 * 
 *
 * @author tim.jeske
 */
public class VCFFilterNodeModel extends NodeModel {
    
	private static final NodeLogger LOGGER = NodeLogger.getLogger(VCFFilterNodeModel.class);
	
	//input
	static final String CFGKEY_VCFIN = "vcf_infile";
	private final SettingsModelString m_vcfin = new SettingsModelString(CFGKEY_VCFIN,"-");
	
	static final String CFGKEY_VEP_SCRIPT = "vepscript";
    private final SettingsModelString m_vepscript = new SettingsModelString(CFGKEY_VEP_SCRIPT,"-");
    
    static final String CFGKEY_VCFTOOLS = "vcf_tools";
    private final SettingsModelString m_vcf_tools = new SettingsModelString(CFGKEY_VCFTOOLS,"-");
	
    //genotype filter
    static final String CFGKEY_FILTER_BY_DP = "filter_by_DP";
    final SettingsModelBoolean m_filter_by_DP = new SettingsModelBoolean(CFGKEY_FILTER_BY_DP,false);
    
    static final String CFGKEY_DP_THRESHOLD = "DP_threshold";
    final SettingsModelInteger m_DP_threshold = new SettingsModelInteger(CFGKEY_DP_THRESHOLD,8);
    
    static final String CFGKEY_FILTER_BY_GQ = "filter_by_GQ";
    final SettingsModelBoolean m_filter_by_GQ = new SettingsModelBoolean(CFGKEY_FILTER_BY_GQ,false);
    
    static final String CFGKEY_GQ_THRESHOLD = "GQ_threshold";
    final SettingsModelInteger m_GQ_threshold = new SettingsModelInteger(CFGKEY_GQ_THRESHOLD,20);
    
    //variant filter
    static final String CFGKEY_FILTER_BY_GQ_MEAN = "filter_by_GQ_MEAN";
    final SettingsModelBoolean m_filter_by_GQ_MEAN = new SettingsModelBoolean(CFGKEY_FILTER_BY_GQ_MEAN,false);
    
    static final String CFGKEY_GQ_MEAN_THRESHOLD = "GQ_mean_threshold";
    final SettingsModelInteger m_GQ_mean_threshold = new SettingsModelInteger(CFGKEY_GQ_MEAN_THRESHOLD,35);
    
    //annotation filter
    static final String CFGKEY_FILTER_ANNOTATIONS = "filter_annotations";
	final SettingsModelBoolean m_filter_annos = new SettingsModelBoolean(CFGKEY_FILTER_ANNOTATIONS,false);
    
	static final String CFGKEY_ANNOTATION="annotation";
    static final String[] ANNOTATIONS_AVAILABLE={"VAT", "VEP"};
    private final SettingsModelString m_annotation = new SettingsModelString(CFGKEY_ANNOTATION, "");
    
    static final String [] SO_TERMS = {"splice_acceptor_variant", "splice_donor_variant",
    	"stop_gained", "frameshift_variant", "stop_lost",
    	"initiator_codon_variant", "inframe_insertion","inframe_deletion", "missense_variant"};
    
    static final String [] DEFAULT_TERMS = {SO_TERMS[0], SO_TERMS[1], SO_TERMS[2], SO_TERMS[3]};
    
    static final String CFGKEY_SO_TERM = "so_term"; 
    
    static final String CFGKEY_TERM_LIST = "term_list";
    
    private final HashSet<String> TERMS	= new HashSet<String>();

    static final String CFGKEY_FILTER = "filter";
    private final SettingsModelString m_filter = new SettingsModelString(CFGKEY_FILTER,"");
    
    
	//output col names
	public static final String OUT_COL1 = "Path2FilteredVCF";
	
	public boolean optionalPort=false;
	
    protected VCFFilterNodeModel() {
    	super(OptionalPorts.createOPOs(1, 1), OptionalPorts.createOPOs(1));
        for(String t: DEFAULT_TERMS) {
        	this.TERMS.add(t);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	//input file
    	String infile;
    	
    	if(optionalPort){	//Input Table available
    		//Get File via Table
    		infile = inData[0].iterator().next().getCell(0).toString();
    	}else{
    		//Get File via FileSelector
    		infile = m_vcfin.getStringValue();
    		if(infile.equals("") || Files.notExists(Paths.get(infile))) {
    			LOGGER.error("No input vcf file specified!");
    		}
    	}
    	
    	String outfile = "";
    	if(m_filter_by_DP.getBooleanValue() || m_filter_by_GQ.getBooleanValue()) {
    		outfile = filterGenotypes(infile, exec);
    		infile = outfile;
    	}
    	if(m_filter_annos.getBooleanValue()) {
    		outfile = filterAnnotations(infile, exec);
    	}
    	
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(outfile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	return new BufferedDataTable[]{outTable};
    }

    private String filterGenotypes(String infile, ExecutionContext exec) throws Exception {
    	
    	String outfile;
    	if(infile.endsWith(".gz")) {
    		outfile = infile.replace("vcf.gz", "genotype.filtered");
    	} else {
    		outfile = infile.replace("vcf", "genotype.filtered");
    	}
    	
    	if(m_filter_by_GQ.getBooleanValue() || m_filter_by_DP.getBooleanValue()) {
    		ArrayList<String> cmd = new ArrayList<>();
    		String vcftools = m_vcf_tools.getStringValue();
    		if(!vcftools.endsWith(System.getProperty("file.separator"))) {
    			vcftools += System.getProperty("file.separator");
    		}
    		vcftools+="vcftools";
    		cmd.add(vcftools);
    		//TODO go on here
    		
    	}
    	
    	return outfile;
    }
    
    private String filterAnnotations(String infile, ExecutionContext exec) throws Exception {
    	String outfile;
    	if(infile.endsWith(".gz")) {
    		outfile = infile.replace("vcf.gz", "annotation.filtered.vcf");
    	} else {
    		outfile = infile.replace("vcf", "annotation.filtered.vcf");
    	}
    	
    	LOGGER.debug("CHOSEN TERMS: "+TERMS);
    	LOGGER.debug("FILTER: "+m_filter.getStringValue());
    	
    	String annotation = m_annotation.getStringValue();
    	LOGGER.debug(annotation);
    	if(annotation.equals("VAT")) {
    		LOGGER.debug("Filter VAT annotation");
    		new VATFilter().filter(infile, outfile, TERMS);
    	} else if(annotation.equals("VEP")) {
    		LOGGER.info("Prepare command for filter_vep.pl");
    		ArrayList<String> cmd = new ArrayList<>();
    		cmd.add("perl");
    		cmd.add(m_vepscript.getStringValue());
    		cmd.add("-i");
    		cmd.add(infile);
    		cmd.add("-o");
    		cmd.add(outfile);
    		
    		String filter_terms = "";
    		Object [] terms =  TERMS.toArray();
    		if(terms.length>0) {
    			cmd.add("--filter");
    			filter_terms += "Consequence is "+terms[0];
    			for(int i=1; i<terms.length; i++) {
        			filter_terms += " or Consequence is "+terms[i];
        		}
    			cmd.add(filter_terms);
    		}
    		
    		String [] filters = m_filter.getStringValue().split(",");
    		if(filters.length>0) {
    			for(String f: filters) {
    				if(f.equals("")) continue;
    				cmd.add("--filter");
    				f = f.trim();
    				cmd.add(f);
    			}
    		}
    		
    		cmd.add("--only_matched");

    		String [] cmd_array = new String[cmd.size()];
    		for(int i = 0; i < cmd.size(); i++) {
    			cmd_array[i] = cmd.get(i);
    		}
    		
    		Executor.executeCommand(cmd_array, exec, LOGGER);
    	}
    	return outfile;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        this.TERMS.clear();
        for(String t: DEFAULT_TERMS) {
        	this.TERMS.add(t);
        }
        optionalPort = false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	String infile_warning = CheckUtils.checkSourceFile(m_vcfin.getStringValue());
    	if(infile_warning != null) {
    		setWarningMessage(infile_warning);
    	}
    	
    	String vep_warning = CheckUtils.checkSourceFile(m_vepscript.getStringValue());
    	if(vep_warning != null) {
    		setWarningMessage(vep_warning);
    	}
    	
    	String vcftools_warning = CheckUtils.checkSourceFile(m_vcf_tools.getStringValue());
    	if(vcftools_warning != null) {
    		setWarningMessage(vcftools_warning);
    	}
    	
    	try{
			inSpecs[0].getColumnNames();
			optionalPort=true;
			
		}catch(NullPointerException e){}
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
	protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_filter_by_DP.saveSettingsTo(settings);
    	m_DP_threshold.saveSettingsTo(settings);
    	m_filter_by_GQ.saveSettingsTo(settings);
    	m_GQ_threshold.saveSettingsTo(settings);
    	m_filter_by_GQ_MEAN.saveSettingsTo(settings);
    	m_GQ_mean_threshold.saveSettingsTo(settings);
    	m_vcf_tools.saveSettingsTo(settings);
		m_filter_annos.saveSettingsTo(settings);
		m_vcfin.saveSettingsTo(settings);
		m_filter.saveSettingsTo(settings);
		m_vepscript.saveSettingsTo(settings);
		m_annotation.saveSettingsTo(settings);
		settings.addStringArray(VCFFilterNodeModel.CFGKEY_TERM_LIST,
				TERMS.toArray(new String[TERMS.size()]));
	}

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_filter_by_DP.loadSettingsFrom(settings);
    	m_DP_threshold.loadSettingsFrom(settings);
    	m_filter_by_GQ.loadSettingsFrom(settings);
    	m_GQ_threshold.loadSettingsFrom(settings);
    	m_filter_by_GQ_MEAN.loadSettingsFrom(settings);
    	m_GQ_mean_threshold.loadSettingsFrom(settings);
    	m_vcf_tools.loadSettingsFrom(settings);
    	m_filter_annos.loadSettingsFrom(settings);
        m_vcfin.loadSettingsFrom(settings);
        m_filter.loadSettingsFrom(settings);
        m_vepscript.loadSettingsFrom(settings);
        m_annotation.loadSettingsFrom(settings);
        
        // clean the old data
    	TERMS.clear();
    	// check, if data is set
        if (settings.containsKey(VCFFilterNodeModel.CFGKEY_TERM_LIST)) {
        	try {
        		// add the values
				for(String s : settings.getStringArray(VCFFilterNodeModel.CFGKEY_TERM_LIST))
					this.TERMS.add(s);
			} catch (InvalidSettingsException e) {
				LOGGER.error(e.getStackTrace());
			}
        } 
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_filter_by_DP.validateSettings(settings);
    	m_DP_threshold.validateSettings(settings);
    	m_filter_by_GQ.validateSettings(settings);
    	m_GQ_threshold.validateSettings(settings);
    	m_filter_by_GQ_MEAN.validateSettings(settings);
    	m_GQ_mean_threshold.validateSettings(settings);
    	m_vcf_tools.validateSettings(settings);
    	m_filter_annos.validateSettings(settings);
        m_vcfin.validateSettings(settings);
        m_filter.validateSettings(settings);
        m_vepscript.validateSettings(settings);
        m_annotation.validateSettings(settings);
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


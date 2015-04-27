package de.helmholtz_muenchen.ibis.ngs.vepfilter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
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
import org.knime.core.node.defaultnodesettings.SettingsModelString;

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
public class LOFFilterNodeModel extends NodeModel {
    
	//input
	static final String CFGKEY_VCFIN = "vcf_infile";
	private final SettingsModelString m_vcfin = new SettingsModelString(CFGKEY_VCFIN,"-");
	
	static final String CFGKEY_ANNOTATION="annotation";
    static final String[] ANNOTATIONS_AVAILABLE={"VAT", "VEP"};
    private final SettingsModelString m_annotation = new SettingsModelString(CFGKEY_ANNOTATION, "");
    
    static final String CFGKEY_VEP_SCRIPT = "vepscript";
    private final SettingsModelString m_vepscript = new SettingsModelString(CFGKEY_VEP_SCRIPT,"-");
	
	//LoF definition
    static final String [] SO_TERMS = {"splice_acceptor_variant", "splice_donor_variant",
    	"stop_gained", "frameshift_variant", "stop_lost",
    	"initiator_codon_variant", "inframe_insertion","inframe_deletion", "missense_variant"};
    
    static final String [] DEFAULT_TERMS = {SO_TERMS[0], SO_TERMS[1], SO_TERMS[2], SO_TERMS[3]};
    
    static final String CFGKEY_SO_TERM = "so_term"; 
    
    static final String CFGKEY_TERM_LIST = "term_list";
    
    private final HashSet<String> TERMS	= new HashSet<String>();
    
    private static final NodeLogger LOGGER = NodeLogger.getLogger(LOFFilterNodeModel.class);
    
    //further filter
    static final String CFGKEY_FILTER = "filter";
    private final SettingsModelString m_filter = new SettingsModelString(CFGKEY_FILTER,"-");
    
	//output col names
	public static final String OUT_COL1 = "Path2FilteredVCF";
	
    protected LOFFilterNodeModel() {
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

    	String infile = m_vcfin.getStringValue();
    	String outfile = infile.replace("vcf", "filtered.vcf");
    	
    	LOGGER.debug("CHOSEN TERMS: "+TERMS);
    	LOGGER.debug("FILTER: "+m_filter.getStringValue());
    	
    	String annotation = m_annotation.getStringValue();
    	LOGGER.debug(annotation);
    	if(annotation.equals("VAT")) {
    		LOGGER.debug("ANNOTATION: VAT");
    	} else if(annotation.equals("VEP")) {
    		System.out.println("here");
    		String cmd = "perl";
    		
    		String vepscript = m_vepscript.getStringValue();
    		if(vepscript.equals("") || Files.notExists(Paths.get(vepscript))) {
    			LOGGER.error("Path to filter_vep.pl not specified!");
    		} else {
    			cmd += " " + vepscript;
    		}
    		
    		if(infile.equals("") || Files.notExists(Paths.get(infile))) {
    			LOGGER.error("Input file does not exist");
    		} else {
    			cmd += " -i " + infile;
    		}
    		
    		cmd += " -o " + outfile;
    		
    		String [] filters = m_filter.getStringValue().split(",");
    		for(String f: filters) {
    			cmd += " --filter \""+f+"\"";
    		}
    		
    		cmd += " --only_matched";
    		System.out.println("hherww");
    		LOGGER.info("execute: "+cmd);
    		Executor.executeCommand(new String[]{cmd}, exec, LOGGER);
    		
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

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        this.TERMS.clear();

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_vcfin.saveSettingsTo(settings);
         m_filter.saveSettingsTo(settings);
         m_vepscript.saveSettingsTo(settings);
         m_annotation.saveSettingsTo(settings);
         settings.addStringArray(LOFFilterNodeModel.CFGKEY_TERM_LIST, TERMS.toArray(new String[TERMS.size()]));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_vcfin.loadSettingsFrom(settings);
        m_filter.loadSettingsFrom(settings);
        m_vepscript.loadSettingsFrom(settings);
        m_annotation.loadSettingsFrom(settings);
        
        // clean the old data
    	TERMS.clear();
    	// check, if data is set
        if (settings.containsKey(LOFFilterNodeModel.CFGKEY_TERM_LIST)) {
        	try {
        		// add the values
				for(String s : settings.getStringArray(LOFFilterNodeModel.CFGKEY_TERM_LIST))
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


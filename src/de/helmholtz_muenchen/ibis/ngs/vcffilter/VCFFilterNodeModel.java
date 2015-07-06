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
public class VCFFilterNodeModel extends NodeModel {
    
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
    
    private static final NodeLogger LOGGER = NodeLogger.getLogger(VCFFilterNodeModel.class);
    
    //further filter
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
    	
    	String outfile;
    	if(infile.endsWith(".gz")) {
    		outfile = infile.replace("vcf.gz", "filtered.vcf");
    	} else {
    		outfile = infile.replace("vcf", "filtered.vcf");
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
//    		String cmd = "perl";
    		ArrayList<String> cmd = new ArrayList<>();
    		cmd.add("perl");
    		
    		String vepscript = m_vepscript.getStringValue();
    		if(vepscript.equals("") || Files.notExists(Paths.get(vepscript))) {
    			LOGGER.error("Path to filter_vep.pl not specified!");
    		} else {
//    			cmd += " " + vepscript;
    			cmd.add(vepscript);
    		}
    		
    		if(infile.equals("") || Files.notExists(Paths.get(infile))) {
    			LOGGER.error("Input file does not exist");
    		} else {
//    			cmd += " -i " + infile;
    			cmd.add("-i");
    			cmd.add(infile);
    		}
    		
//    		cmd += " -o " + outfile;
    		cmd.add("-o");
    		cmd.add(outfile);
    		
//    		String filter_terms = " --filter \"";
    		
    		String filter_terms = "";
    		Object [] terms =  TERMS.toArray();
    		if(terms.length>0) {
    			cmd.add("--filter");
    			filter_terms += "Consequence is "+terms[0];
    			for(int i=1; i<terms.length; i++) {
        			filter_terms += " or Consequence is "+terms[i];
        		}
    			cmd.add(filter_terms);
//    			cmd += filter_terms;
    		}
    		
    		String [] filters = m_filter.getStringValue().split(",");
    		if(filters.length>0) {
    			for(String f: filters) {
    				if(f.equals("")) continue;
    				cmd.add("--filter");
    				f = f.trim();
    				cmd.add(f);
//    				cmd += " --filter "+f;
    			}
    		}
    		
//    		cmd += " --only_matched";
    		cmd.add("--only_matched");

    		String [] cmd_array = new String[cmd.size()];
    		for(int i = 0; i < cmd.size(); i++) {
    			cmd_array[i] = cmd.get(i);
    		}
    		
//    		Executor.executeCommand(new String[]{cmd}, exec, LOGGER);
    		Executor.executeCommand(cmd_array, exec, LOGGER);
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
         m_vcfin.saveSettingsTo(settings);
         m_filter.saveSettingsTo(settings);
         m_vepscript.saveSettingsTo(settings);
         m_annotation.saveSettingsTo(settings);
         settings.addStringArray(VCFFilterNodeModel.CFGKEY_TERM_LIST, TERMS.toArray(new String[TERMS.size()]));
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


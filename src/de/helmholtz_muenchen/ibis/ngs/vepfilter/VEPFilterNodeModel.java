package de.helmholtz_muenchen.ibis.ngs.vepfilter;

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
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;
/**
 * This is the model implementation of VEPFilter.
 * 
 *
 * @author Tim Jeske
 */
public class VEPFilterNodeModel extends HTExecutorNodeModel {
    
	//configuration keys
    static final String CFGKEY_VEP_SCRIPT = "vepscript";
    static final String CFGKEY_SO_TERM = "so_term"; 
    static final String CFGKEY_TERM_LIST = "term_list";
    static final String CFGKEY_FILTER = "filter";    
    static final String CFGKEY_OUTFOLDER = "outfolder";
	static final String CFGKEY_OVERWRITE = "overwrite";
	
	//settings models
	private final SettingsModelString m_vepscript = new SettingsModelString(CFGKEY_VEP_SCRIPT,"-");
    private final SettingsModelOptionalString m_filter = new SettingsModelOptionalString(CFGKEY_FILTER,"",false);
	private final SettingsModelString m_outfolder = new SettingsModelString(CFGKEY_OUTFOLDER,"");
	private final SettingsModelBoolean m_overwrite = new SettingsModelBoolean(CFGKEY_OVERWRITE, false);
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(VEPFilterNodeModel.class);
	
	public static final String OUT_COL1 = "Path2FilteredVCF";
	
	static final String [] SO_TERMS = {"splice_acceptor_variant", "splice_donor_variant",
	    	"stop_gained", "frameshift_variant", "stop_lost",
	    	"initiator_codon_variant", "inframe_insertion","inframe_deletion", "missense_variant"};
	    
	static final String [] DEFAULT_TERMS = {SO_TERMS[0], SO_TERMS[1], SO_TERMS[2], SO_TERMS[3]};
	
	private int vcf_index;
	private final HashSet<String> TERMS	= new HashSet<String>();
	
    protected VEPFilterNodeModel() {
    	super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
        for(String t: DEFAULT_TERMS) {
        	this.TERMS.add(t);
        }
        m_vepscript.setStringValue(IBISKNIMENodesPlugin.getDefault().getToolPathPreference("filter_vep.pl"));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String infile = inData[0].iterator().next().getCell(vcf_index).toString();
    	
    	if(Files.notExists(Paths.get(infile))) {
    		throw new InvalidSettingsException("Input VCF does not exist!");
    	}
    	
    	String outfile = filterAnnotations(infile, exec);
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(outfile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	return new BufferedDataTable[]{outTable};
    }
    
    private String filterAnnotations(String infile, ExecutionContext exec) throws Exception {
    	
		String outfile;
		if (infile.endsWith(".gz")) {
			outfile = infile.replace(".vcf.gz", ".VEPfiltered.vcf");
		} else {
			outfile = IO.replaceFileExtension(infile,
					".VEPfiltered.vcf");
		}
		
		String outfolder = m_outfolder.getStringValue();
    	if(outfolder!=null && !outfolder.equals("")) {
    		outfile = outfolder + System.getProperty("file.separator") + new File(outfile).getName();
    	}

		LOGGER.debug("CHOSEN TERMS: " + TERMS);
		LOGGER.debug("FILTER: " + m_filter.getStringValue());

		LOGGER.info("Prepare command for filter_vep.pl");
		ArrayList<String> cmd = new ArrayList<>();
		cmd.add("perl");
		cmd.add(m_vepscript.getStringValue());
		cmd.add("-i");
		cmd.add(infile);
		cmd.add("-o");
		cmd.add(outfile);

		String filter_terms = "";
		Object[] terms = TERMS.toArray();
		if (terms.length > 0) {
			cmd.add("--filter");
			filter_terms += "Consequence is " + terms[0];
			for (int i = 1; i < terms.length; i++) {
				filter_terms += " or Consequence is " + terms[i];
			}
			cmd.add(filter_terms);
		}

		String[] filters = m_filter.getStringValue().split(",");
		if (filters.length > 0) {
			for (String f : filters) {
				if (f.equals(""))
					continue;
				cmd.add("--filter");
				f = f.trim();
				cmd.add(f);
			}
		}

		cmd.add("--only_matched");
		
		if(m_overwrite.getBooleanValue()) {
			cmd.add("--force_overwrite");
		}

		String[] cmd_array = new String[cmd.size()];
		for (int i = 0; i < cmd.size(); i++) {
			cmd_array[i] = cmd.get(i);
		}

		File lockFile = new File(outfile + SuccessfulRunChecker.LOCK_ENDING);
		super.executeCommand(cmd_array, exec, lockFile);

		return outfile;
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
    	
    	vcf_index=-1;
    	
    	String vep_warning = CheckUtils.checkSourceFile(m_vepscript.getStringValue());
    	if(vep_warning != null) {
    		setWarningMessage(vep_warning);
    	}
    	
    	for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
    		if(inSpecs[0].getColumnSpec(i).getType().toString().equals("VCFCell")) {
    			vcf_index = i;
    		}
    	}
    	
    	if(vcf_index==-1) {
    		throw new InvalidSettingsException("This node is not compatible with the precedent node as there is no VCF file in the output table!");
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
		m_filter.saveSettingsTo(settings);
		m_vepscript.saveSettingsTo(settings);
		m_outfolder.saveSettingsTo(settings);
		m_overwrite.saveSettingsTo(settings);
		settings.addStringArray(VEPFilterNodeModel.CFGKEY_TERM_LIST,
				TERMS.toArray(new String[TERMS.size()]));
	}

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
        m_filter.loadSettingsFrom(settings);
        m_vepscript.loadSettingsFrom(settings);
        m_outfolder.loadSettingsFrom(settings);
        m_overwrite.loadSettingsFrom(settings);
        
    	TERMS.clear();
        if (settings.containsKey(VEPFilterNodeModel.CFGKEY_TERM_LIST)) {
        	try {
				for(String s : settings.getStringArray(VEPFilterNodeModel.CFGKEY_TERM_LIST))
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
    	super.validateSettings(settings);
        m_filter.validateSettings(settings);
        m_vepscript.validateSettings(settings);
        m_outfolder.validateSettings(settings);
        m_overwrite.validateSettings(settings);
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


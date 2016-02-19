package de.helmholtz_muenchen.ibis.ngs.vcftoolsfilter;

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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;
/**
 * This is the model implementation of LOFFilter.
 * 
 *
 * @author Tim Jeske
 */
public class VCFToolsFilterNodeModel extends HTExecutorNodeModel {
    
	private static final NodeLogger LOGGER = NodeLogger.getLogger(VCFToolsFilterNodeModel.class);
	
    
    static final String CFGKEY_VCFTOOLS = "vcf_tools";
    private final SettingsModelString m_vcf_tools = new SettingsModelString(CFGKEY_VCFTOOLS,"-");
    
    static final String CFGKEY_FILL_AN_AC = "fill_an_ac";
    final SettingsModelBoolean m_fill_an_ac = new SettingsModelBoolean(CFGKEY_FILL_AN_AC,false);
    
    //genotype filter
    static final String CFGKEY_FILTER_BY_DP = "filter_by_DP";
    final SettingsModelBoolean m_filter_by_DP = new SettingsModelBoolean(CFGKEY_FILTER_BY_DP,false);
    
    static final String CFGKEY_DP_THRESHOLD = "DP_threshold";
    final SettingsModelInteger m_DP_threshold = new SettingsModelInteger(CFGKEY_DP_THRESHOLD,8);
    
    static final String CFGKEY_FILTER_BY_AD = "filter_by_AD";
    final SettingsModelBoolean m_filter_by_AD = new SettingsModelBoolean(CFGKEY_FILTER_BY_AD,false);
    
    static final String CFGKEY_AD_THRESHOLD = "AD_threshold";
    final SettingsModelInteger m_AD_threshold = new SettingsModelInteger(CFGKEY_AD_THRESHOLD,8);
    
    static final String CFGKEY_FILTER_BY_GQ = "filter_by_GQ";
    final SettingsModelBoolean m_filter_by_GQ = new SettingsModelBoolean(CFGKEY_FILTER_BY_GQ,false);
    
    static final String CFGKEY_GQ_THRESHOLD = "GQ_threshold";
    final SettingsModelIntegerBounded m_GQ_threshold = new SettingsModelIntegerBounded(CFGKEY_GQ_THRESHOLD,20,0,99);
    
    //variant filter
    static final String CFGKEY_FILTER_PASS = "filter_pass";
    final SettingsModelBoolean m_filter_pass = new SettingsModelBoolean(CFGKEY_FILTER_PASS,false);
    
    static final String CFGKEY_FILTER_CALL_RATE = "filter_callRate";
    final SettingsModelBoolean m_filter_by_callRate = new SettingsModelBoolean(CFGKEY_FILTER_CALL_RATE,false);
    
    static final String CFGKEY_CALL_RATE_THRESHOLD = "callRate_threshold";
    final SettingsModelDoubleBounded m_callRate_threshold = new SettingsModelDoubleBounded(CFGKEY_CALL_RATE_THRESHOLD, 0.88,0.0,1.0);
    
    static final String CFGKEY_FILTER_BY_GQ_MEAN = "filter_by_GQ_MEAN";
    final SettingsModelBoolean m_filter_by_GQ_MEAN = new SettingsModelBoolean(CFGKEY_FILTER_BY_GQ_MEAN,false);
    
    static final String CFGKEY_GQ_MEAN_THRESHOLD = "GQ_mean_threshold";
    final SettingsModelIntegerBounded m_GQ_mean_threshold = new SettingsModelIntegerBounded(CFGKEY_GQ_MEAN_THRESHOLD,35,0,99);
    
    
	//output col names
	public static final String OUT_COL1 = "Path2FilteredVCF";
	
	private int vcf_index;
	
    protected VCFToolsFilterNodeModel() {
    	super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {


    	String infile = inData[0].iterator().next().getCell(vcf_index).toString();
    	
    	if(Files.notExists(Paths.get(infile))) {
    		throw new InvalidSettingsException("Input VCF file does not exist!");
    	}
    	
    	String outfile = infile;
    	if(m_filter_by_DP.getBooleanValue() || m_filter_by_GQ.getBooleanValue() || m_filter_by_AD.getBooleanValue()) {
    		outfile = filterGenotypes(infile, exec);
    		infile = outfile;
    	}
    	if(m_filter_pass.getBooleanValue() || m_filter_by_GQ_MEAN.getBooleanValue() || m_filter_by_callRate.getBooleanValue()) {
    		outfile = filterVariants(infile, exec);
    		infile = outfile;
    	}
    	if(m_fill_an_ac.getBooleanValue()) {
    		outfile = postprocess(infile,exec);
    	}
    	
    	
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

    private String filterGenotypes(String infile, ExecutionContext exec) throws Exception {
    	
    	boolean is_gz = infile.endsWith(".gz");
    	String outfile;
    	if(is_gz) {
    		outfile = infile.replace(".vcf.gz", ".genotype_filtered");
    	} else {
    		outfile = IO.replaceFileExtension(infile, ".genotype_filtered");
    	}
    	
    	if(m_filter_by_GQ.getBooleanValue() || m_filter_by_DP.getBooleanValue()) {
    		ArrayList<String> cmd = new ArrayList<>();
    		String vcftools = m_vcf_tools.getStringValue();
    		cmd.add(vcftools);
    		if(is_gz) {
    			cmd.add("--gzvcf");
    		} else {
    			cmd.add("--vcf");
    		}
    		cmd.add(infile);
    		cmd.add("--out");
    		cmd.add(outfile);
    		cmd.add("--minDP");
    		cmd.add(m_DP_threshold.getIntValue()+"");
    		cmd.add("--minGQ");
    		cmd.add(m_GQ_threshold.getIntValue()+"");
    		cmd.add("--recode");
    		cmd.add("--recode-INFO-all");
    		
    		String [] cmd_array = new String[cmd.size()];
    		for(int i = 0; i < cmd.size(); i++) {
    			cmd_array[i] = cmd.get(i);
    		}
    		
    		outfile = outfile+".recode.vcf";
    		File lockFile = new File(outfile+SuccessfulRunChecker.LOCK_ENDING);
    		super.executeCommand(cmd_array, exec, lockFile);
    	}
    	
    	if(m_filter_by_AD.getBooleanValue()) {
    		
    		if(infile.endsWith(".gz")) {
    			LOGGER.warn("AD filter can only process unzipped vcf files!");
    			return outfile;
    		}
    		
    		outfile = infile.replace(".vcf",".AD_filtered.vcf");
    		
    		ArrayList<String> cmd = new ArrayList<>();
			cmd.add("perl");
			cmd.add(IO.getScriptPath()+"scripts/perl/filterByAD.pl");
			cmd.add(infile);
			cmd.add(m_AD_threshold.getIntValue()+"");
			cmd.add("./.");
			
			String[] cmd_array = new String[cmd.size()];
			for (int i = 0; i < cmd.size(); i++) {
				cmd_array[i] = cmd.get(i);
			}

			File lockFile = new File(outfile + SuccessfulRunChecker.LOCK_ENDING);

			super.executeCommand(cmd_array, exec, lockFile);
    		
    	}
    	
    	return outfile;
    }
    
    private String filterVariants(String infile, ExecutionContext exec) throws Exception {
    	String outfile = "";
    	
    	if(m_filter_pass.getBooleanValue() || m_filter_by_callRate.getBooleanValue()) {
    		boolean is_gz = infile.endsWith(".gz");
    		if(is_gz) {
        		outfile = infile.replace(".vcf.gz", ".variant_filtered");
        	} else {
        		outfile = IO.replaceFileExtension(infile, ".variant_filtered");
        	}
    		
    		ArrayList<String> cmd = new ArrayList<>();
    		String vcftools = m_vcf_tools.getStringValue();
    		cmd.add(vcftools);
    		if(is_gz) {
    			cmd.add("--gzvcf");
    		} else {
    			cmd.add("--vcf");
    		}
    		cmd.add(infile);
    		cmd.add("--out");
    		cmd.add(outfile);
    		if(m_filter_pass.getBooleanValue()) {
    			cmd.add("--remove-filtered-all");
    		}
    		if(m_filter_by_callRate.getBooleanValue()) {
    			cmd.add("--max-missing");
    			cmd.add(m_callRate_threshold.getDoubleValue()+"");
    		}
    		
    		cmd.add("--recode");
    		cmd.add("--recode-INFO-all");
    		
    		String [] cmd_array = new String[cmd.size()];
    		for(int i = 0; i < cmd.size(); i++) {
    			cmd_array[i] = cmd.get(i);
    		}
    		
    		outfile = outfile+".recode.vcf";
    		File lockFile = new File(outfile+SuccessfulRunChecker.LOCK_ENDING);
    		
    		super.executeCommand(cmd_array, exec, lockFile);
    		infile = outfile;
    	}
    	
    	if(m_filter_by_GQ_MEAN.getBooleanValue()) {
    		
    		if(infile.endsWith(".gz")) {
    			LOGGER.warn("Mean GQ filter can only process unzipped vcf files!");
    			return outfile;
    		}
    		
    		outfile = infile.replace(".vcf",".GQFiltered.vcf");
    		
    		ArrayList<String> cmd = new ArrayList<>();
			cmd.add("perl");
			cmd.add(IO.getScriptPath()+"scripts/perl/filterGQMean.pl");
			cmd.add(infile);
			cmd.add(m_GQ_mean_threshold.getIntValue()+"");
			cmd.add("0");
			
			String[] cmd_array = new String[cmd.size()];
			for (int i = 0; i < cmd.size(); i++) {
				cmd_array[i] = cmd.get(i);
			}

			File lockFile = new File(outfile + SuccessfulRunChecker.LOCK_ENDING);

			super.executeCommand(cmd_array, exec, lockFile);
    	}
    	return outfile;
    }
    
    private String postprocess(String infile, ExecutionContext exec) throws Exception {
    	boolean is_gz = infile.endsWith(".gz");
    	String outfile;
    	if(is_gz) {
    		outfile = infile.replace(".vcf.gz", ".adjusted_AN_AC.vcf");
    	} else {
    		outfile = infile.replace(".vcf", ".adjusted_AN_AC.vcf");
    	}
    	
    	String cmd = "";
    	if(is_gz) {
			cmd = "zcat";
		} else {
			cmd = "cat";
		}
		
		cmd += " " + infile;
		
    	if(m_fill_an_ac.getBooleanValue()) {
    		
    		String vcftools = m_vcf_tools.getStringValue();
    		int pos = vcftools.lastIndexOf(System.getProperty("file.separator"));
    		vcftools = vcftools.substring(0,pos)+ System.getProperty("file.separator")+"fill-an-ac";

    		cmd += " | " + vcftools;
    	}

    	File lockFile = new File(outfile+SuccessfulRunChecker.LOCK_ENDING);	
    	cmd += " > " + outfile;
    	super.executeCommand(new String[]{"/usr/bin/bash","-c",cmd},exec,lockFile);
    	
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
    	
    	boolean a = m_filter_by_DP.getBooleanValue();
    	boolean b = m_filter_by_GQ.getBooleanValue();
    	boolean c = m_filter_by_callRate.getBooleanValue();
    	boolean d = m_fill_an_ac.getBooleanValue();
    	boolean f = m_filter_pass.getBooleanValue();
    	vcf_index = -1;
    	
    	if(a || b || c || d || f) {
    		String vcftools_warning = CheckUtils.checkSourceFile(m_vcf_tools.getStringValue());
        	if(vcftools_warning != null) {
        		setWarningMessage(vcftools_warning);
        	}
    	}
    	
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
    	m_filter_by_DP.saveSettingsTo(settings);
    	m_filter_by_AD.saveSettingsTo(settings);
    	m_DP_threshold.saveSettingsTo(settings);
    	m_filter_by_GQ.saveSettingsTo(settings);
    	m_GQ_threshold.saveSettingsTo(settings);
    	m_filter_pass.saveSettingsTo(settings);
    	m_filter_by_callRate.saveSettingsTo(settings);
    	m_callRate_threshold.saveSettingsTo(settings);
    	m_filter_by_GQ_MEAN.saveSettingsTo(settings);
    	m_GQ_mean_threshold.saveSettingsTo(settings);
    	m_vcf_tools.saveSettingsTo(settings);
		m_fill_an_ac.saveSettingsTo(settings);
	}

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
    	m_filter_by_DP.loadSettingsFrom(settings);
    	m_filter_by_AD.loadSettingsFrom(settings);
    	m_DP_threshold.loadSettingsFrom(settings);
    	m_filter_by_GQ.loadSettingsFrom(settings);
    	m_GQ_threshold.loadSettingsFrom(settings);
    	m_filter_pass.loadSettingsFrom(settings);
    	m_filter_by_callRate.loadSettingsFrom(settings);
    	m_callRate_threshold.loadSettingsFrom(settings);
    	m_filter_by_GQ_MEAN.loadSettingsFrom(settings);
    	m_GQ_mean_threshold.loadSettingsFrom(settings);
    	m_vcf_tools.loadSettingsFrom(settings);
        m_fill_an_ac.loadSettingsFrom(settings); 
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.validateSettings(settings);
    	m_filter_by_DP.validateSettings(settings);
    	m_filter_by_AD.validateSettings(settings);
    	m_DP_threshold.validateSettings(settings);
    	m_filter_by_GQ.validateSettings(settings);
    	m_GQ_threshold.validateSettings(settings);
    	m_filter_by_callRate.validateSettings(settings);
    	m_callRate_threshold.validateSettings(settings);
    	m_filter_pass.validateSettings(settings);
    	m_filter_by_GQ_MEAN.validateSettings(settings);
    	m_GQ_mean_threshold.validateSettings(settings);
    	m_vcf_tools.validateSettings(settings);
        m_fill_an_ac.validateSettings(settings);
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


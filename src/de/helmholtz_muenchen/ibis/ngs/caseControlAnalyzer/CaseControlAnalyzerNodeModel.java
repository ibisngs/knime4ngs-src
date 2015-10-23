package de.helmholtz_muenchen.ibis.ngs.caseControlAnalyzer;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.TreeMap;

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
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of CaseControlAnalyzer.
 * 
 *
 * @author Tim Jeske
 */
public class CaseControlAnalyzerNodeModel extends NodeModel {
    
	//model file configure keys
	static final String CFGKEY_MODEL_GENE_ID = "model_gene_id";
	static final String CFGKEY_FREQ = "freq";
	static final String CFGKEY_PSEUDO_FREQ = "pseudo_freq";
	
	//summary file configure keys
	static final String CFGKEY_SUMMARY_GENE_ID = "summary_gene_id";
	static final String CFGKEY_CASE_COND = "case_cond";
	static final String CFGKEY_CASE_NCOND = "case_ncond";
	static final String CFGKEY_CONTROL_COND = "control_cond";
	static final String CFGKEY_CONTROL_NCOND = "control_ncond";
	
	//model file settings models
    private final SettingsModelString m_model_gene_id = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_MODEL_GENE_ID, "gene_id");
    private final SettingsModelString m_freq = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_FREQ, "lof_freq");
    private final SettingsModelDouble m_pseudo_freq = new SettingsModelDouble(CaseControlAnalyzerNodeModel.CFGKEY_PSEUDO_FREQ,0.0);
    
    //summary file settings models
    private final SettingsModelString m_summary_gene_id = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_SUMMARY_GENE_ID, "gene_id");
    private final SettingsModelString m_case_cond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_COND, "aff_case");
    private final SettingsModelString m_case_ncond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_NCOND, "un_case");
    private final SettingsModelString m_control_cond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CONTROL_COND, "aff_ctrl");
    private final SettingsModelString m_control_ncond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_NCOND, "un_ctrl");

    private static final String OUT_COL1 = "Path2ExtendedSummary";
    
    private boolean model_given = false;
    
    protected CaseControlAnalyzerNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(OptionalPorts.createOPOs(2, 2), OptionalPorts.createOPOs(1));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String summary_file, model_file, outfile;
    	HashMap<String, Double> gene2frequency;
    	HashMap<String, ContingencyTable> gene2table;
    	LinkedHashMap<String, HashMap<String, Double>> gene2pvalues;
    	LinkedList<String> ordered_genes;
    	
    	summary_file = inData[0].iterator().next().getCell(0).toString();
    	gene2table = readSummaryFile(summary_file, m_summary_gene_id.getStringValue(), m_case_cond.getStringValue(), m_case_ncond.getStringValue(), m_control_cond.getStringValue(), m_control_ncond.getStringValue());
    	
    	gene2frequency = new HashMap<>();
    	if(model_given) {
    		model_file = inData[1].iterator().next().getCell(0).toString();
    		gene2frequency = readModelFile(model_file, m_model_gene_id.getStringValue(), m_freq.getStringValue());
    	}
    	
    	gene2pvalues = new LinkedHashMap<>();
    	
    	HashMap<String, Double> fisher = getFisherPvalues(gene2table);
    	gene2pvalues.put("fisher_exact", fisher);
    	gene2pvalues.put("fisher_exact_adj", adjustPvalues(fisher));
    	//if order by fisher
    	ordered_genes = getOrder(fisher);
    	
    	HashMap<String, Double> bg = getBinomialBackgroundPvalues(gene2table, gene2frequency, m_pseudo_freq.getDoubleValue());
    	gene2pvalues.put("binomial_background", bg);
    	gene2pvalues.put("fisher_exact_adj", adjustPvalues(bg));
    	//if order by fisher
    	ordered_genes = getOrder(bg);
    	
    	outfile= IO.replaceFileExtension(summary_file, "extended.tsv");
    	
    	writeResults(outfile, summary_file, gene2pvalues,ordered_genes);
    	
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
   
    private HashMap<String, Double> readModelFile(String file, String gene_id_header, String freq_header) {
    	HashMap<String, Double> result = new HashMap<>();
    	//TODO implement
    	return result;
    }
    
    private HashMap<String, ContingencyTable> readSummaryFile(String file, String gene_id, String case_cond, String case_ncond, String control_cond, String control_ncond) {
    	HashMap<String, ContingencyTable> result = new HashMap<>();
    	//TODO implement
    	return result;
    }
    
    private HashMap<String, Double> getFisherPvalues(HashMap<String,ContingencyTable> gene2contingency) {
    	HashMap<String, Double> result = new HashMap<>();
    	//TODO implement
    	return result;
    }
    
    private HashMap<String, Double> getBinomialBackgroundPvalues(HashMap<String,ContingencyTable> gene2contingency, HashMap<String, Double> gene2freq, double pseudo_freq) {
    	HashMap<String, Double> result = new HashMap<>();
    	//TODO implement
    	return result;
    }
    
    private HashMap<String, Double> adjustPvalues(HashMap<String, Double> pvalues) {
    	HashMap<String, Double> result = new HashMap<>();
    	//TODO implement
    	return result;
    }
    
    private LinkedList<String> getOrder(HashMap<String, Double> pvalues) {
    	LinkedList<String> genes = new LinkedList<>();
    	
    	ValueComparator vc = new ValueComparator(pvalues);

		TreeMap<String,Double> sorted_gene_statistic = new TreeMap<String, Double>(vc);

		sorted_gene_statistic.putAll(pvalues);
		genes.addAll(sorted_gene_statistic.keySet());
    	
    	return genes;
    }
    
	class ValueComparator implements Comparator<String> {

	    HashMap<String, Double> base;
	    public ValueComparator(HashMap<String, Double> base) {
	        this.base = base;
	    }
   
	    public int compare(String a, String b) {
	        return Double.compare(base.get(a),base.get(b));
	    }
	}
     
    private void  writeResults(String outfile, String summary_file, LinkedHashMap<String, HashMap<String,Double>> pvalues, LinkedList<String> order) {
    	
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

    	model_given = false;
    	
    	if(inSpecs[1] != null) {
    		model_given = true;
    		checkSettingsModelString(m_model_gene_id);
        	checkSettingsModelString(m_freq);
    	}
    	
    	m_model_gene_id.setEnabled(model_given);
    	m_freq.setEnabled(model_given);
    	
    	checkSettingsModelString(m_summary_gene_id);
    	checkSettingsModelString(m_case_cond);
    	checkSettingsModelString(m_case_ncond);
    	checkSettingsModelString(m_control_cond);
    	checkSettingsModelString(m_control_ncond);
    	
    	Double tmp = m_pseudo_freq.getDoubleValue();
    	if(tmp==null || tmp <0.0) {
    		throw new InvalidSettingsException("The value for the pseudo frequency must be greater than or equal to 0.0");
    	}
    	
        return new DataTableSpec[]{new DataTableSpec(
        		new DataColumnSpec[]{
        				new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()})};
    }
    
    private void checkSettingsModelString (SettingsModelString m) throws InvalidSettingsException {
    	String tmp = m.getStringValue();
    	if(tmp.equals("") || tmp==null) {
    		throw new InvalidSettingsException("The "+ m.getKey()+" field must not be empty!");
    	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_model_gene_id.saveSettingsTo(settings);
         m_freq.saveSettingsTo(settings);
         m_pseudo_freq.saveSettingsTo(settings);
         m_summary_gene_id.saveSettingsTo(settings);
         m_case_cond.saveSettingsTo(settings);
         m_case_ncond.saveSettingsTo(settings);
         m_control_cond.saveSettingsTo(settings);
         m_control_ncond.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_model_gene_id.loadSettingsFrom(settings);
    	m_freq.loadSettingsFrom(settings);
    	m_pseudo_freq.loadSettingsFrom(settings);
    	m_summary_gene_id.loadSettingsFrom(settings);
    	m_case_cond.loadSettingsFrom(settings);
    	m_case_ncond.loadSettingsFrom(settings);
    	m_control_cond.loadSettingsFrom(settings);
    	m_control_ncond.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_model_gene_id.validateSettings(settings);
        m_freq.validateSettings(settings);
        m_pseudo_freq.validateSettings(settings);
        m_summary_gene_id.validateSettings(settings);
        m_case_cond.validateSettings(settings);
        m_case_ncond.validateSettings(settings);
        m_control_cond.validateSettings(settings);
        m_control_ncond.validateSettings(settings);
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


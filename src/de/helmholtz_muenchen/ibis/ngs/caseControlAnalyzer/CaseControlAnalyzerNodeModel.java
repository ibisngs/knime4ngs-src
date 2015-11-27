package de.helmholtz_muenchen.ibis.ngs.caseControlAnalyzer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
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
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
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
	
	//summary file configure keys
	static final String CFGKEY_SUMMARY_GENE_ID = "summary_gene_id";
	static final String CFGKEY_CASE_COND = "case_cond";
	static final String CFGKEY_CASE_NCOND = "case_ncond";
	static final String CFGKEY_CONTROL_COND = "control_cond";
	static final String CFGKEY_CONTROL_NCOND = "control_ncond";
	static final String CFGKEY_IGNORE = "ignore";
	
	//method configure keys
	static final String CFGKEY_FISHER = "fisher_exact";
	static final String CFGKEY_BINOMIAL_BACKGROUND = "binomial_background";
	static final String CFGKEY_PSEUDO_FREQ = "pseudo_freq";
	static final String CFGKEY_HYPER_BACKGROUND = "hypergeometric_background";
	static final String CFGKEY_POP_SIZE = "pop_size";
	static final String CFGKEY_ORDER_BY = "order_by";
	static final String [] METHODS = {CFGKEY_FISHER, CFGKEY_BINOMIAL_BACKGROUND, CFGKEY_HYPER_BACKGROUND};
	
	//model file settings models
    private final SettingsModelString m_model_gene_id = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_MODEL_GENE_ID, "gene_id");
    private final SettingsModelString m_freq = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_FREQ, "lof_freq");
    
    //summary file settings models
    private final SettingsModelString m_summary_gene_id = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_SUMMARY_GENE_ID, "gene_id");
    private final SettingsModelString m_case_cond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_COND, "aff_case");
    private final SettingsModelString m_case_ncond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CASE_NCOND, "un_case");
    private final SettingsModelString m_control_cond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CONTROL_COND, "aff_ctrl");
    private final SettingsModelString m_control_ncond = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_CONTROL_NCOND, "un_ctrl");
    private final SettingsModelBoolean m_ignore = new SettingsModelBoolean(CaseControlAnalyzerNodeModel.CFGKEY_IGNORE, true);
    
    //method settings models
    private final SettingsModelBoolean m_fisher = new SettingsModelBoolean(CaseControlAnalyzerNodeModel.CFGKEY_FISHER,true);
    private final SettingsModelBoolean m_bin_back = new SettingsModelBoolean(CaseControlAnalyzerNodeModel.CFGKEY_BINOMIAL_BACKGROUND,true);
    private final SettingsModelDoubleBounded m_pseudo_freq = new SettingsModelDoubleBounded(CaseControlAnalyzerNodeModel.CFGKEY_PSEUDO_FREQ,0.0,0.0,1.0);
    private final SettingsModelBoolean m_hyper = new SettingsModelBoolean(CaseControlAnalyzerNodeModel.CFGKEY_HYPER_BACKGROUND,true);
    private final SettingsModelIntegerBounded m_pop_size = new SettingsModelIntegerBounded(CaseControlAnalyzerNodeModel.CFGKEY_POP_SIZE, 1, 1, Integer.MAX_VALUE);
    private final SettingsModelString m_order_by = new SettingsModelString(CaseControlAnalyzerNodeModel.CFGKEY_ORDER_BY,METHODS[0]);
    
    private static final String OUT_COL1 = "Path2ExtendedSummary";
    
    private boolean model_given = false;
    
	protected static final NodeLogger logger = NodeLogger.getLogger(CaseControlAnalyzerNodeModel.class);

    
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

    	Statistics stats = new Statistics();
    	
    	String summary_file, model_file, outfile;
    	HashMap<String, Double> gene2frequency;
    	HashMap<String, ContingencyTable> gene2table;
    	LinkedHashMap<String, HashMap<String, Double>> gene2pvalues;
    	LinkedList<String> ordered_genes = new LinkedList<>();
    	
    	summary_file = inData[0].iterator().next().getCell(0).toString();
    	gene2table = readSummaryFile(summary_file, m_summary_gene_id.getStringValue(), m_case_cond.getStringValue(), m_case_ncond.getStringValue(), m_control_cond.getStringValue(), m_control_ncond.getStringValue());
    	
    	gene2frequency = new HashMap<>();
    	if(model_given) {
    		model_file = inData[1].iterator().next().getCell(0).toString();
    		gene2frequency = readModelFile(model_file, m_model_gene_id.getStringValue(), m_freq.getStringValue());
    	}
    	
    	gene2pvalues = new LinkedHashMap<>();
    	
    	logger.debug("Compute Fisher statistic");
    	if(m_fisher.getBooleanValue()) {
        	HashMap<String, Double> fisher = getFisherPvalues(gene2table,stats);
    		gene2pvalues.put("fisher_exact", fisher);
    		gene2pvalues.put("fisher_exact_adj", adjustPvaluesFDR(fisher,stats));
    		if(m_order_by.getStringValue().equals(METHODS[0])) {
    			ordered_genes = getOrder(fisher);
    		}
    	}
    	
    	logger.debug("Compute binomial background statistic");
    	if(m_bin_back.getBooleanValue() && model_given) {
        	HashMap<String, Double> bg = getBinomialBackgroundPvalues(gene2table, gene2frequency, m_pseudo_freq.getDoubleValue(), stats);
    		gene2pvalues.put("binomial_background", bg);
    		gene2pvalues.put("binomial_background_adj", adjustPvaluesFDR(bg, stats));
    		if(m_order_by.getStringValue().equals(METHODS[1])) {
    			ordered_genes = getOrder(bg);
    		}
    	}
    	
    	logger.debug("Compute hypergeometric background statistic");
    	if(m_hyper.getBooleanValue() && model_given) {
    		HashMap<String, Double> hg = getHypergeometricBackgroundPvalues(gene2table, gene2frequency, m_pop_size.getIntValue(), stats);
    		gene2pvalues.put("hypergeometric_background", hg);
    		gene2pvalues.put("hypergeometric_background_adj", adjustPvaluesFDR(hg, stats));
    		if(m_order_by.getStringValue().equals(METHODS[2])) {
    			ordered_genes = getOrder(hg);
    		}
    	}
    	
//    	gene2pvalues.put("normal_approx", getNormalApprox(gene2table, gene2frequency,m_pop_size.getIntValue(),stats));
    	
    	outfile= IO.replaceFileExtension(summary_file, "extended.tsv");
    	
    	//order by first map in gene2pvalues if ordering was incorrectly specified
    	if(ordered_genes.size()==0) {
    		ordered_genes = getOrder(gene2pvalues.get(gene2pvalues.keySet().iterator().next()));
    	}
    	
    	stats.quit();
    	writeResults(outfile, summary_file, m_summary_gene_id.getStringValue(), gene2pvalues,gene2frequency,ordered_genes);
    	
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
   


	private HashMap<String, Double> readModelFile(String file, String gene_id_header, String freq_header) throws IOException, InvalidSettingsException {
    	HashMap<String, Double> result = new HashMap<>();
    	
    	int gene_index = -1;
    	int freq_index = -1;
    	
    	String line;
    	String [] fields;
    	BufferedReader br = Files.newBufferedReader(Paths.get(file));
    	
    	String header = br.readLine();
    	if(header.contains(gene_id_header) && header.contains(freq_header)) {
    		String [] cols = header.split("\t");
    		for(int i = 0; i < cols.length; i++) {
    			if(cols[i].equals(gene_id_header)) {
    				gene_index = i;
    			} else if (cols[i].equals(freq_header)) {
    				freq_index = i;
    			}
    		}
    	} else {
    		throw new InvalidSettingsException("Specified column headers are not found in model file!");
    	}
    	
    	while((line=br.readLine())!=null) {
    		fields = line.split("\t");
    		result.put(fields[gene_index], Double.parseDouble(fields[freq_index]));
    	}
    	
    	br.close();
    	
    	return result;
    }
    
    private HashMap<String, ContingencyTable> readSummaryFile(String file, String gene_id, String case_cond, String case_ncond, String control_cond, String control_ncond) throws IOException, InvalidSettingsException {
    	HashMap<String, ContingencyTable> result = new HashMap<>();
    	
    	int gene_index = -1;
    	int a_index = -1 ;
    	int b_index = -1;
    	int c_index = -1;
    	int d_index = -1;
    	
    	String line;
    	String [] fields;
    	BufferedReader br = Files.newBufferedReader(Paths.get(file));
    	
    	String header = br.readLine();
    	if(header.contains(gene_id) && header.contains(case_cond) && header.contains(case_ncond) && header.contains(control_cond) && header.contains(control_ncond)) {
    		String [] cols = header.split("\t");
    		for(int i = 0; i < cols.length; i++) {
    			if(cols[i].equals(gene_id)) {
    				gene_index = i;
    			} else if(cols[i].equals(case_cond)) {
    				a_index = i;
    			} else if(cols[i].equals(case_ncond)) {
    				b_index = i;
    			} else if(cols[i].equals(control_cond)) {
    				c_index = i;
    			} else if(cols[i].equals(control_ncond)) {
    				d_index = i;
    			}
    		}
    	} else {
    		throw new InvalidSettingsException("Specified column headers are not found in summary file!");
    	}
    	
    	int case_aff, case_un, control_aff, control_un;
    	double frac_case, frac_control;
    	
    	while((line=br.readLine())!=null) {
    		fields = line.split("\t");
    		case_aff = Integer.parseInt(fields[a_index]);
    		case_un = Integer.parseInt(fields[b_index]);
    		control_aff = Integer.parseInt(fields[c_index]);
    		control_un = Integer.parseInt(fields[d_index]);
    		
    		if(m_ignore.getBooleanValue()) {
    			frac_case = (double)case_aff/(double)(case_aff+case_un);
    			frac_control = (double)control_aff/(double)(control_aff+control_un);
    			if(Double.compare(frac_case, frac_control) <= 0) {
    				continue;
    			}
    		}
    		
    		result.put(fields[gene_index], new ContingencyTable(case_aff, case_un, control_aff, control_un));
    	}
    	
    	br.close();
    	return result;
    }
    
    private HashMap<String, Double> getFisherPvalues(HashMap<String,ContingencyTable> gene2contingency, Statistics stats) {
    	HashMap<String, Double> result = new HashMap<>();
    	for(String gene: gene2contingency.keySet()) {
    		result.put(gene, stats.getFisherOneTailedGreater(gene2contingency.get(gene)));
    	}
    	return result;
    }
    
    private HashMap<String, Double> getBinomialBackgroundPvalues(HashMap<String,ContingencyTable> gene2contingency, HashMap<String, Double> gene2freq, double pseudo_freq, Statistics stats) {
    	HashMap<String, Double> result = new HashMap<>();
		double freq;
    	for(String gene: gene2contingency.keySet()) {
			
    		if(gene2freq.containsKey(gene)) {
    			freq = (gene2freq.get(gene) +pseudo_freq)/(1.0+pseudo_freq);
    		} else {
    			freq = pseudo_freq;
    		}

    		result.put(gene, stats.getBinomialBackground(gene2contingency.get(gene), freq));
		}
    	return result;
    }
    
    private HashMap<String, Double> getHypergeometricBackgroundPvalues(HashMap<String,ContingencyTable> gene2contingency, HashMap<String, Double> gene2freq, int pop_size, Statistics stats) {
    	HashMap<String, Double> result = new HashMap<>();
		int pop_cases = 0;
    	for(String gene: gene2contingency.keySet()) {
    		pop_cases = 0;
    		if(gene2freq.containsKey(gene)) {
    			pop_cases = (int)Math.round(gene2freq.get(gene)*pop_size);
    		}
    		result.put(gene, stats.getHypergeometricBackground(gene2contingency.get(gene), pop_size, pop_cases));
		}
    	return result;
    }
    
//    private HashMap<String, Double> getNormalApprox(HashMap<String, ContingencyTable> gene2table,
//			HashMap<String, Double> gene2frequency, int pop_size, Statistics stats) {
//		HashMap<String, Double> result = new HashMap<>();
//		int pop_cases = 0;
//		for(String gene: gene2table.keySet()) {
//    		pop_cases = 0;
//    		if(gene2frequency.containsKey(gene)) {
//    			pop_cases = (int)Math.round(gene2frequency.get(gene)*pop_size);
//    		}
//    		result.put(gene, stats.getNormalApproximation(gene2table.get(gene), pop_size, pop_cases));
//		}
//		return result;
//	}
    
    private HashMap<String, Double> adjustPvaluesFDR(HashMap<String, Double> pvalues, Statistics stats) {
    	HashMap<String, Double> result = new HashMap<>();

    	LinkedList<String> genes = getOrder(pvalues);
    	
    	double [] a = new double[genes.size()];
    	
    	for(int i = 0; i<genes.size(); i++) {
    		a[i] = pvalues.get(genes.get(i));
    	}
    	
    	a = stats.adjustP(a, "fdr");
    	
    	for(int i = 0; i<genes.size(); i++) {
    		result.put(genes.get(i), a[i]);
    	}

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
	    	if(base.get(a).isNaN() || base.get(a) > base.get(b)) {
	    		return 1;
	    	}
	        return -1;
	    }
	}
     
    private void  writeResults(String outfile, String summary_file, String gene_id_header, LinkedHashMap<String, HashMap<String,Double>> pvalues, HashMap<String,Double> gene2frequency, LinkedList<String> order) throws IOException {
    	HashMap<String, String> content = new HashMap<>();
    	BufferedReader br = Files.newBufferedReader(Paths.get(summary_file));
    	
    	//read header
    	String header = br.readLine() + "\tbackground_freq";
    	int gene_index = -1;
    	String [] cols = header.split("\t");
    	for(int i = 0; i < cols.length; i++) {
    		if(cols[i].equals(gene_id_header)) {
    			gene_index = i;
    		}
    	}
    	
    	for(String stat : pvalues.keySet()) {
    		header += "\t" + stat;
    	}
    	
    	//read input file and save in content map
    	String line;
    	while((line=br.readLine())!= null) {
    		content.put(line.split("\t")[gene_index], line);
    	}
    	br.close();
    	
    	BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
    	bw.write(header);
    	bw.newLine();

    	for(String gene: order) {
    		line = content.get(gene);
    		if(gene2frequency.containsKey(gene)) {
    			line +="\t"+(gene2frequency.get(gene)+m_pseudo_freq.getDoubleValue())/(1.0 + m_pseudo_freq.getDoubleValue());
    		} else {
    			line += "\t"+m_pseudo_freq.getDoubleValue();
    		}
    		for(HashMap<String,Double> map: pvalues.values()) {
    			line += "\t" + map.get(gene);
    		}
    		bw.write(line);
    		bw.newLine();
    	}
    	
    	bw.close();
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
         m_fisher.saveSettingsTo(settings);
         m_bin_back.saveSettingsTo(settings);
         m_hyper.saveSettingsTo(settings);
         m_pop_size.saveSettingsTo(settings);
         m_order_by.saveSettingsTo(settings);
         m_ignore.saveSettingsTo(settings);
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
    	m_fisher.loadSettingsFrom(settings);
    	m_bin_back.loadSettingsFrom(settings);
    	m_hyper.loadSettingsFrom(settings);
    	m_pop_size.loadSettingsFrom(settings);
    	m_order_by.loadSettingsFrom(settings);
    	m_ignore.loadSettingsFrom(settings);
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
        m_fisher.validateSettings(settings);
        m_bin_back.validateSettings(settings);
        m_hyper.validateSettings(settings);
        m_pop_size.loadSettingsFrom(settings);
        m_order_by.validateSettings(settings);
        m_ignore.validateSettings(settings);
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


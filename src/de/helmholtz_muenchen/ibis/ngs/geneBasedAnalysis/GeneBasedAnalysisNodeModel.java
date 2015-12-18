package de.helmholtz_muenchen.ibis.ngs.geneBasedAnalysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;

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
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer.CaseControlAnalyzerNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.ContingencyTable;
import de.helmholtz_muenchen.ibis.utils.ngs.Statistics;

/**
 * This is the model implementation of CaseControlAnalyzer.
 * 
 *
 * @author Tim Jeske
 */
public class GeneBasedAnalysisNodeModel extends CaseControlAnalyzerNodeModel {
	
	//method configure keys
	static final String CFGKEY_FISHER = "fisher_exact";
	static final String CFGKEY_BINOMIAL_BACKGROUND = "binomial_background";
	static final String CFGKEY_PSEUDO_FREQ = "pseudo_freq";
	static final String CFGKEY_HYPER_BACKGROUND = "hypergeometric_background";
	static final String CFGKEY_POP_SIZE = "pop_size";
	static final String CFGKEY_ORDER_BY = "order_by";
	static final String [] METHODS = {CFGKEY_FISHER, CFGKEY_BINOMIAL_BACKGROUND, CFGKEY_HYPER_BACKGROUND};
    
    //method settings models
    private final SettingsModelBoolean m_fisher = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_FISHER,true);
    private final SettingsModelBoolean m_bin_back = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_BINOMIAL_BACKGROUND,true);
    private final SettingsModelDoubleBounded m_pseudo_freq = new SettingsModelDoubleBounded(GeneBasedAnalysisNodeModel.CFGKEY_PSEUDO_FREQ,0.0,0.0,1.0);
    private final SettingsModelBoolean m_hyper = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_HYPER_BACKGROUND,true);
    private final SettingsModelIntegerBounded m_pop_size = new SettingsModelIntegerBounded(GeneBasedAnalysisNodeModel.CFGKEY_POP_SIZE, 1, 1, Integer.MAX_VALUE);
    private final SettingsModelString m_order_by = new SettingsModelString(GeneBasedAnalysisNodeModel.CFGKEY_ORDER_BY,METHODS[0]);
    
	protected static final NodeLogger logger = NodeLogger.getLogger(CaseControlAnalyzerNodeModel.class);

	private String outfile;
	
    protected GeneBasedAnalysisNodeModel() {
        super();
    }
    
	@Override
	protected void performAnalysis(BufferedDataTable[] inData, ExecutionContext exec,
			HashMap<String, Double> gene2frequency, int pop_size, HashMap<String, ContingencyTable> gene2table, String gene_id) throws IOException {
		Statistics stats = new Statistics();
    	
    	String summary_file;
    	
    	String [] genes;
    	ContingencyTable [] tables;
    	double [] frequencies = null;
    	double [] fisher = null;
    	double [] fisher_adj = null;
    	double [] binomial = null;
    	double [] binomial_adj = null;
    	double [] hyper = null;
    	double [] hyper_adj = null;
    	
    	summary_file = inData[0].iterator().next().getCell(0).toString();
    	outfile = IO.replaceFileExtension(summary_file, ".extended.tsv");
    	
    	//create arrays
    	int n = gene2table.size();
    	genes = new String[n];
    	tables = new ContingencyTable[n];
    	frequencies = new double[n];
    	
    	Object[] g = gene2table.keySet().toArray();
    	String my_gene;
    	for(int i = 0; i < n; i++) {
    		my_gene = (String)g[i];
    		genes[i] = my_gene;
    		tables[i] = gene2table.get(my_gene);
    		if(gene2frequency.containsKey(my_gene)) {
    			frequencies[i] = gene2frequency.get(my_gene);
    		} else {
    			frequencies[i] = 0.0;
    		}
    	}
    	
    	logger.debug("Compute Fisher statistic");
    	if(m_fisher.getBooleanValue()) {
        	fisher = stats.getFisherOneTailedGreater(tables);
        	fisher_adj = stats.adjustP(fisher, "fdr");
    	}
    	
    	logger.debug("Compute binomial background statistic");
    	if(m_bin_back.getBooleanValue()) {
    		binomial = stats.getBinomialBackground(tables, frequencies, m_pseudo_freq.getDoubleValue());
        	binomial_adj = stats.adjustP(binomial, "fdr");
    	}
    	
    	logger.debug("Compute hypergeometric background statistic");
    	if(m_hyper.getBooleanValue()) {
    		hyper = stats.getHypergeometricBackground(tables, pop_size, frequencies);
    		hyper_adj = stats.adjustP(hyper, "fdr");
    	}
    	
    	stats.quit();
    	
    	//prepare output
    	//get order
    	Integer [] indices = new Integer [genes.length];
    	for(int i = 0; i < indices.length; i++) {
    		indices[i] = i;
    	}

    	String order_by = m_order_by.getStringValue();
    	ValueComparator vc = null;
    	if(order_by.equals(METHODS[0]) && fisher != null) {
    		vc = new ValueComparator(fisher);
    	} else if (order_by.equals(METHODS[1]) && binomial != null) {
    		vc = new ValueComparator(binomial);
    	} else if (order_by.equals(METHODS[2]) && hyper != null) {
    		vc = new ValueComparator(hyper);
    	}
    	
    	Arrays.sort(indices, vc);
    	
    	//get header and arrays to print
    	ArrayList<double[]> p_values = new ArrayList<>();
    	ArrayList<String> headers = new ArrayList<>();
    	
    	if(fisher!=null) {
    		p_values.add(fisher);
    		headers.add("fisher");
    		p_values.add(fisher_adj);
    		headers.add("fisher_adj");
    	}
    	
    	if(binomial!=null) {
    		p_values.add(binomial);
    		headers.add("binomial");
    		p_values.add(binomial_adj);
    		headers.add("biomial_adj");
    	}
    	
    	if(hyper!=null) {
    		p_values.add(hyper);
    		headers.add("hyper");
    		p_values.add(hyper_adj);
    		headers.add("hyper_adj");
    	}
    	
    	outfile= IO.replaceFileExtension(summary_file, "extended.tsv");
    	writeResults(outfile, summary_file, gene_id, headers, p_values, indices, genes, frequencies);
	}

	@Override
	protected String getOutCol() {
		return "Path2ExtendedSummary";
	}

	@Override
	protected String getOutfile() {
		return this.outfile;
	}
    
	class ValueComparator implements Comparator<Integer> {

	    private final double [] array;
	    public ValueComparator(double [] arr) {
	        this.array = arr; 
	    }

	    @Override
	    public int compare(Integer index1, Integer index2) {
	         if(new Double(array[index1]).isNaN() || array[index1] > array[index2]) {
	        	 return 1; 
	         }
	         return -1;
	    }
	}
    
     
    private void  writeResults(String outfile, String summary_file, String gene_id_header, ArrayList<String> headers, ArrayList<double[]> p_values, Integer[] indices, String [] genes, double [] frequencies) throws IOException {
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
    	
    	for(String stat : headers) {
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

    	for(int i :indices) {
    		line = content.get(genes[i]);
    		line += "\t"+frequencies[i];
    		for(double [] stat: p_values) {
    			line += "\t" + stat[i];
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
	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
		 m_pseudo_freq.saveSettingsTo(settings);
         m_fisher.saveSettingsTo(settings);
         m_bin_back.saveSettingsTo(settings);
         m_hyper.saveSettingsTo(settings);
         m_pop_size.saveSettingsTo(settings);
         m_order_by.saveSettingsTo(settings);
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings) throws InvalidSettingsException {
		m_pseudo_freq.loadSettingsFrom(settings);
    	m_fisher.loadSettingsFrom(settings);
    	m_bin_back.loadSettingsFrom(settings);
    	m_hyper.loadSettingsFrom(settings);
    	m_pop_size.loadSettingsFrom(settings);
    	m_order_by.loadSettingsFrom(settings);
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings) throws InvalidSettingsException {
		m_pseudo_freq.validateSettings(settings);
        m_fisher.validateSettings(settings);
        m_bin_back.validateSettings(settings);
        m_hyper.validateSettings(settings);
        m_pop_size.loadSettingsFrom(settings);
        m_order_by.validateSettings(settings);
	}

}
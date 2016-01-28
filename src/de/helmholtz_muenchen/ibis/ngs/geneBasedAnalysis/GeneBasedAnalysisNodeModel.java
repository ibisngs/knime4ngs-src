package de.helmholtz_muenchen.ibis.ngs.geneBasedAnalysis;

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
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer.*;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer.Identifier.*;
import de.helmholtz_muenchen.ibis.utils.ngs.ContingencyTable;
import de.helmholtz_muenchen.ibis.utils.ngs.Statistics;

/**
 * This is the model implementation of CaseControlAnalyzer.
 * 
 *
 * @author Tim Jeske
 */
public class GeneBasedAnalysisNodeModel extends CaseControlAnalyzerNodeModel {
	
	//arrays for string selection fields
	static final String [] RESOLUTION = {"gene","transcript"};
	static final String [] ALTERNATIVES = {"two.sided", "greater", "less"};
	
	//method configure keys
	static final String CFGKEY_RESOLUTION = "resolution";

	static final String CFGKEY_FISHER = "fisher_exact";
	static final String CFGKEY_ALT_FISHER = "fisher_alternative";
	static final String CFGKEY_WILCOXON = "wilcoxon";
	static final String CFGKEY_ALT_WILCOXON = "wilcoxon_alternative";
//	static final String CFGKEY_BINOMIAL_BACKGROUND = "binomial_background";
//	static final String CFGKEY_PSEUDO_FREQ = "pseudo_freq";
	static final String CFGKEY_HYPER_BACKGROUND = "hypergeometric_background";
	static final String CFGKEY_POP_SIZE = "pop_size";
	static final String CFGKEY_ORDER_BY = "order_by";
	static final String [] METHODS = {CFGKEY_FISHER, CFGKEY_WILCOXON, CFGKEY_HYPER_BACKGROUND};
	
    
    //method settings models
	private final SettingsModelString m_resolution = new SettingsModelString(GeneBasedAnalysisNodeModel.CFGKEY_RESOLUTION, GeneBasedAnalysisNodeModel.RESOLUTION[0]);
    private final SettingsModelBoolean m_fisher = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_FISHER,true);
    private final SettingsModelString m_alt_fisher = new SettingsModelString(GeneBasedAnalysisNodeModel.CFGKEY_ALT_FISHER,"two.sided");
    private final SettingsModelBoolean m_wilcoxon = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_WILCOXON,true);
    private final SettingsModelString m_alt_wilcoxon = new SettingsModelString(GeneBasedAnalysisNodeModel.CFGKEY_ALT_WILCOXON,"two.sided");
//    private final SettingsModelBoolean m_bin_back = new SettingsModelBoolean(GeneBasedAnalysisNodeModel.CFGKEY_BINOMIAL_BACKGROUND,true);
//    private final SettingsModelDoubleBounded m_pseudo_freq = new SettingsModelDoubleBounded(GeneBasedAnalysisNodeModel.CFGKEY_PSEUDO_FREQ,0.0,0.0,1.0);
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
			HashMap<String, Double> gene2frequency, int pop_size, MatrixSummary ms) throws IOException {
		
		HashMap<String, ContingencyTable> entity2table;
		HashMap<String, CaseControlArray> entity2array;
		if(m_resolution.getStringValue().equals("gene")) {
			entity2table = ms.toTables(new GeneIdentifier());
			entity2array = ms.toArrays(new GeneIdentifier());
		} else {
			entity2table = ms.toTables();
			entity2array = ms.toArrays();
		}
		
		Statistics stats = new Statistics();
    	
    	String [] genes;
    	ContingencyTable [] tables;
    	CaseControlArray [] arrays;
    	double [] frequencies = null;
    	double [] fisher = null;
    	double [] fisher_adj = null;
    	double [] wilcoxon = null;
    	double [] wilcoxon_adj = null;
//    	double [] binomial = null;
//    	double [] binomial_adj = null;
    	double [] hyper = null;
    	double [] hyper_adj = null;
    	
    	String summary_file = inData[0].iterator().next().getCell(0).toString();
    	outfile = IO.replaceFileExtension(summary_file, "."+m_resolution.getStringValue()+"_analysis.tsv");
    	
    	//create arrays
    	int n = entity2table.size();
    	genes = new String[n];
    	tables = new ContingencyTable[n];
    	arrays = new CaseControlArray[n];
    	frequencies = new double[n];
    	
    	Object[] g = entity2table.keySet().toArray();
    	String entity;
    	for(int i = 0; i < n; i++) {
    		entity = (String)g[i];
    		genes[i] = entity;
    		tables[i] = entity2table.get(entity);
    		arrays[i] = entity2array.get(entity);
    		if(gene2frequency.containsKey(entity.split("_")[0])) { //works in case of genes and transcripts
    			frequencies[i] = gene2frequency.get(entity.split("_")[0]);
    		} else {
    			frequencies[i] = 0.0;
    		}
    	}
    	
    	logger.debug("Compute Fisher statistic");
    	if(m_fisher.getBooleanValue()) {
        	fisher = stats.getFisherTest(tables,m_alt_fisher.getStringValue());
        	fisher_adj = stats.adjustP(fisher, "fdr");
    	}
    	
    	logger.debug("Compute Wilcoxon statistic");
    	if(m_wilcoxon.getBooleanValue()) {
        	wilcoxon = stats.getWilcoxon(arrays,m_alt_wilcoxon.getStringValue());
        	wilcoxon_adj = stats.adjustP(wilcoxon, "fdr");
    	}
    	
    	
//    	logger.debug("Compute binomial background statistic");
//    	if(m_bin_back.getBooleanValue()) {
//    		binomial = stats.getBinomialBackground(tables, frequencies, m_pseudo_freq.getDoubleValue());
//        	binomial_adj = stats.adjustP(binomial, "fdr");
//    	}
    	
    	logger.debug("Compute hypergeometric background statistic");
    	if(m_hyper.getBooleanValue()) {
    		hyper = stats.getHypergeometricBackground(tables, pop_size, frequencies);
    		hyper_adj = stats.adjustP(hyper, "fdr");
    	}
    	
    	stats.quit();
    	
    	
    	//get header and arrays to print
    	ArrayList<double[]> p_values = new ArrayList<>();
    	ArrayList<String> headers = new ArrayList<>();
    	
    	if(fisher!=null) {
    		p_values.add(fisher);
    		headers.add("fisher");
    		p_values.add(fisher_adj);
    		headers.add("fisher_adj");
    	}
    	
    	if(wilcoxon != null) {
    		p_values.add(wilcoxon);
    		headers.add("wilcoxon");
    		p_values.add(wilcoxon_adj);
    		headers.add("wilcoxon_adj");
    	}
    	
//    	if(binomial!=null) {
//    		p_values.add(binomial);
//    		headers.add("binomial");
//    		p_values.add(binomial_adj);
//    		headers.add("biomial_adj");
//    	}
    	
    	if(hyper!=null) {
    		p_values.add(hyper);
    		headers.add("hyper");
    		p_values.add(hyper_adj);
    		headers.add("hyper_adj");
    	}
    	
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
    	} else if (order_by.equals(METHODS[1]) && wilcoxon != null) {
    		vc = new ValueComparator(wilcoxon);
//    	} else if (order_by.equals(METHODS[2]) && binomial != null) {
//    		vc = new ValueComparator(binomial);
    	} else if (order_by.equals(METHODS[2]) && hyper != null) {
    		vc = new ValueComparator(hyper);
    	}
    	
    	//wrong user input
    	if(vc == null) {
    		vc = new ValueComparator(p_values.get(0));
    	}
    	
    	Arrays.sort(indices, vc);
    	
    	writeResults(outfile, headers, p_values, indices, genes, entity2table, frequencies);
	}

	@Override
	protected String getOutCol() {
		return "Path2GeneSummary";
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
	    	 if (new Double(array[index1]).isNaN() && new Double(array[index2]).isNaN()) {
	        	 return 0;
	         } else if(new Double(array[index1]).isNaN() || array[index1] > array[index2]) {
	        	 return 1; 
	         } 
	         return -1;
	    }
	}
     
    private void  writeResults(String outfile, ArrayList<String> headers, ArrayList<double[]> p_values, Integer[] indices, String [] genes,  HashMap<String, ContingencyTable> gene2table, double [] frequencies) throws IOException {
    	
    	String entity;
    	if(m_resolution.getStringValue().equals("gene")) {
    		entity = "gene_id\tgene_symbol";
    	} else {
    		entity = "transcript\tgene_id\tgene_symbol";
    	}
    	
    	String header = entity + "\taff_case\taff_ctrl\tun_case\tun_ctrl\tbackground_freq";
    	for(String stat : headers) {
    		header += "\t" + stat;
    	}
    	
    	BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
    	bw.write(header);
    	bw.newLine();

    	String line;
    	for(int i :indices) {
    		line = genes[i].replaceAll("_", "\t")+"\t"+gene2table.get(genes[i]).verticalToString()+"\t"+frequencies[i];
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
		m_resolution.saveSettingsTo(settings);
//		m_pseudo_freq.saveSettingsTo(settings);
		m_wilcoxon.saveSettingsTo(settings);
		m_alt_wilcoxon.saveSettingsTo(settings);
        m_fisher.saveSettingsTo(settings);
        m_alt_fisher.saveSettingsTo(settings);
//      m_bin_back.saveSettingsTo(settings);
        m_hyper.saveSettingsTo(settings);
        m_pop_size.saveSettingsTo(settings);
        m_order_by.saveSettingsTo(settings);
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings) throws InvalidSettingsException {
		m_resolution.loadSettingsFrom(settings);
//		m_pseudo_freq.loadSettingsFrom(settings);
		m_wilcoxon.loadSettingsFrom(settings);
		m_alt_wilcoxon.loadSettingsFrom(settings);
    	m_fisher.loadSettingsFrom(settings);
    	m_alt_fisher.loadSettingsFrom(settings);
//    	m_bin_back.loadSettingsFrom(settings);
    	m_hyper.loadSettingsFrom(settings);
    	m_pop_size.loadSettingsFrom(settings);
    	m_order_by.loadSettingsFrom(settings);
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings) throws InvalidSettingsException {
		m_resolution.validateSettings(settings);
//		m_pseudo_freq.validateSettings(settings);
		m_wilcoxon.validateSettings(settings);
		m_alt_wilcoxon.validateSettings(settings);
        m_fisher.validateSettings(settings);
        m_alt_fisher.loadSettingsFrom(settings);
//      m_bin_back.validateSettings(settings);
        m_hyper.validateSettings(settings);
        m_pop_size.loadSettingsFrom(settings);
        m_order_by.validateSettings(settings);
	}

}
package de.helmholtz_muenchen.ibis.ngs.geneSetAnalysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.ngs.BioEntity;
import de.helmholtz_muenchen.ibis.utils.ngs.ContingencyTable;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer.*;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer.Identifier.*;
import de.helmholtz_muenchen.ibis.utils.ngs.Statistics;

/**
 * This is the model implementation of GeneSetAnalysis.
 * 
 *
 * @author Tim Jeske
 */
public class GeneSetAnalysisNodeModel extends CaseControlAnalyzerNodeModel {
	
	static final String CFGKEY_GENE_SET_INFILE = "gene_set_infile";
	static final String CFGKEY_GENE_SET = "gene_set_name";
		
	final SettingsModelString m_genesetin = new SettingsModelString(GeneSetAnalysisNodeModel.CFGKEY_GENE_SET_INFILE,"");
	final SettingsModelString m_genesetname = new SettingsModelString(GeneSetAnalysisNodeModel.CFGKEY_GENE_SET,"gene_set");
	
	protected static final NodeLogger logger = NodeLogger.getLogger(GeneSetAnalysisNodeModel.class);
	
	private String outfile;
	
    protected GeneSetAnalysisNodeModel() {
    	super();
    }

	@Override
	protected void performAnalysis(BufferedDataTable[] inData, ExecutionContext exec,
			HashMap<String, Double> gene2frequency, int pop_size, MatrixSummary ms)
					throws IOException, InvalidSettingsException {
		
    	
    	
    	// check/read gene set file
    	logger.debug("Read gene set file...");
    	String geneset_file = m_genesetin.getStringValue();
    	if(geneset_file.equals("") || Files.notExists(Paths.get(geneset_file))) {
    		throw new InvalidSettingsException("No gene set file specified!");
    	}

    	//generate map: gene_symbol -> ContingencyTable
    	HashMap<String, ContingencyTable> gene2table = ms.toTables(new EntityIdentifier(BioEntity.GENE_SYMBOL));
    	
    	//generate maps: gene_set(s) <-> gene_symbol(s)
    	HashMap<String, HashSet<String>> set2genes = readGeneSetFile(geneset_file, gene2table.keySet());
    	HashMap<String, HashSet<String>> genes2sets = new HashMap<>();
    	HashSet<String> tmp;
    	for(String s: set2genes.keySet()) {
    		for(String g: set2genes.get(s)) {
    			if(genes2sets.containsKey(g)) {
    				genes2sets.get(g).add(s);
    			} else {
    				tmp = new HashSet<>();
    				tmp.add(s);
    				genes2sets.put(g, tmp);
    			}
    		}
    	}
    	
    	//compute Wilcoxon test
    	HashMap<String,Double> set2pval = new HashMap<>();
    	
    	ArrayList<ContingencyTable> my_cts;
    	ContingencyTable [] cts;
    	ArrayList<Double> my_bg_freq;
    	double [] bg_freqs;

    	Statistics stats = new Statistics();
    	
    	for(String s: set2genes.keySet()) {
    		my_cts = new ArrayList<>();
    		my_bg_freq = new ArrayList<>();
    		for(String gene: set2genes.get(s)) {
    			my_cts.add(gene2table.get(gene));
    			if(gene2frequency.containsKey(gene)) {
    				my_bg_freq.add(gene2frequency.get(gene));
    			} else {
    				my_bg_freq.add(0.0);
    			}
    		}
    		cts = new ContingencyTable[my_cts.size()];
    		bg_freqs = new double[my_bg_freq.size()];
    		for(int i = 0; i < my_cts.size(); i++) {
    			cts[i] = my_cts.get(i);
    		}
    		for(int i = 0; i < my_bg_freq.size(); i++) {
    			bg_freqs[i] = my_bg_freq.get(i);
    		}
    		if(bg_freqs.length >= 1 && cts.length >= 1) {
    			set2pval.put(s, stats.getWilcoxonSignedRankTest(cts,pop_size,bg_freqs));
    		}
    	}
    	
    	//sort p-values
    	ValueComparator vc = new ValueComparator(set2pval);
    	TreeMap<String,Double> sorted_pvals = new TreeMap<>(vc);
    	sorted_pvals.putAll(set2pval);
    	
    	//adjust p-values
    	double [] pvals = new double[sorted_pvals.size()];
    	int i = 0;
    	for(Entry<String,Double> e: sorted_pvals.entrySet()) {
    		pvals[i] = e.getValue();
    		i++;
    	}
    	
    	double [] adj_pvals = stats.adjustP(pvals, "fdr");
    	stats.quit();
    	
    	//generate map: gene_set -> ContingencyTable
    	HashMap<String, ContingencyTable> set2table = ms.toTables(new GeneSetIdentifier(genes2sets));
    	
    	
    	String matrix_file = inData[0].iterator().next().getCell(0).toString();
    	
    	String set_name = m_genesetname.getStringValue();
    	if(set_name.equals("")) {
    		set_name = "gene_set";
    	}
    	
    	outfile = IO.replaceFileExtension(matrix_file, "."+set_name+"_analysis.tsv");
    	writeResults(outfile, set2genes, set2table, sorted_pvals, adj_pvals); 
	}

	
	

	
	class ValueComparator implements Comparator<String> {

	    private final HashMap<String,Double> p_vals;
	    public ValueComparator(HashMap<String,Double> p_vals) {
	        this.p_vals = p_vals; 
	    }

	    @Override
	    public int compare(String a, String b) {
	         if(p_vals.get(a).isNaN()|| p_vals.get(a) > p_vals.get(b)) {
	        	 return 1; 
	         }
	         return -1;
	    }
	}
	
	private HashMap<String, HashSet<String>> readGeneSetFile(String file, Set<String> genes) throws IOException {
		HashMap<String, HashSet<String>> set2genes = new HashMap<>();
		
		if(file==null) {
			return set2genes;
		}
		
		BufferedReader br = Files.newBufferedReader(Paths.get(file));
		String line,gene;
		String [] fields;
		
		while((line=br.readLine())!=null) {
			line = line.trim();
			fields = line.split("\t");
			HashSet<String> tmp = new HashSet<>();
			for(int i = 2; i < fields.length; i++) {
				gene = fields[i].toUpperCase();
				if(genes.contains(gene)) {
					tmp.add(gene);
				} else {
					logger.debug("No gene information found for: "+fields[i] +" in set "+fields[0]);
				}
			}
			set2genes.put(fields[0], tmp);
		}
		
		br.close();
		
		return set2genes;
	}
	
    private void writeResults(String outfile, HashMap<String, HashSet<String>> set2genes, HashMap<String, ContingencyTable> set2table, TreeMap<String, Double> set2pval, double[] adj_pvals) throws IOException {

    	String header = "set\tsize\taff_case\taff_ctrl\tun_case\tun_ctrl\tp_Wilcoxon\tp_adj\taff_genes";
    	
    	BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
    	bw.write(header);
    	bw.newLine();
    	
    	int i = 0;
    	for(Entry<String, Double> s: set2pval.entrySet()) {
    		bw.write(s.getKey()+"\t"+set2genes.get(s.getKey()).size()+"\t"+set2table.get(s.getKey()).verticalToString()+"\t"+s.getValue()+"\t"+adj_pvals[i]+"\t"+set2genes.get(s.getKey()));
    		bw.newLine();
    		i++;
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



	@Override
	protected String getOutCol() {
		return "Path2GeneSetSummary";
	}

	@Override
	protected String getOutfile() {
		return this.outfile;
	}

	@Override
	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
		m_genesetin.saveSettingsTo(settings);
		m_genesetname.saveSettingsTo(settings);
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings) throws InvalidSettingsException {
		m_genesetin.loadSettingsFrom(settings);
		m_genesetname.loadSettingsFrom(settings);
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings) throws InvalidSettingsException {
		m_genesetin.validateSettings(settings);
		m_genesetname.validateSettings(settings);
	}
}
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
import java.util.TreeMap;

import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.ngs.ContingencyTable;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.caseControlAnalyzer.CaseControlAnalyzerNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.Statistics;

/**
 * This is the model implementation of GeneSetAnalysis.
 * 
 *
 * @author Tim Jeske
 */
public class GeneSetAnalysisNodeModel extends CaseControlAnalyzerNodeModel {
	
	static final int SET_ID_INDEX = 0;
	static final int SET_INDEX = 6;
		
	protected static final NodeLogger logger = NodeLogger.getLogger(GeneSetAnalysisNodeModel.class);
	
	private String outfile;
	
    protected GeneSetAnalysisNodeModel() {
    	super();
    }

	@Override
	protected void performAnalysis(BufferedDataTable[] inData, ExecutionContext exec,
			HashMap<String, Double> gene2frequency, int pop_size, HashMap<String, ContingencyTable> gene2table, String gene_id)
					throws IOException {
		//get input
    	String gene_set_sum_file = inData[0].iterator().next().getCell(1).toString();
    	outfile = IO.replaceFileExtension(gene_set_sum_file, ".extended.tsv");
    	
    	//store input and results
    	HashMap<String,HashSet<String>> set2genes = readGeneSetSummary(gene_set_sum_file);
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
    		if(bg_freqs.length > 1 && cts.length > 1) {
    			set2pval.put(s, stats.getWilcoxonSignedRankTest(cts,pop_size,bg_freqs));
    		}
    	}
    	
    	writeResults(gene_set_sum_file, outfile, set2pval); 
	}

	class ValueComparator implements Comparator<String> {

	    private final HashMap<String,Double> p_vals;
	    public ValueComparator(HashMap<String,Double> p_vals) {
	        this.p_vals = p_vals; 
	    }

	    @Override
	    public int compare(String a, String b) {
	         if(p_vals.get(a) > p_vals.get(b)) {
	        	 return 1; 
	         }
	         return -1;
	    }
	}
	
	private HashMap<String, HashSet<String>> readGeneSetSummary(String gene_set_sum_file) throws IOException {
		HashMap<String, HashSet<String>> res = new HashMap<>();
    	
    	BufferedReader br = Files.newBufferedReader(Paths.get(gene_set_sum_file));
		String line;
		String [] fields;
		HashSet<String> tmp;
		while((line=br.readLine())!=null) {
			fields = line.split("\t");
			tmp = new HashSet<>();
			for(String s: fields[SET_INDEX].split(",")) {
				tmp.add(s.trim());
			}
			res.put(fields[SET_ID_INDEX],tmp);
		}
		return res;
	}
	
    private void writeResults(String gene_set_sum_file, String outfile, HashMap<String, Double> set2pval) throws IOException {
    	HashMap<String,String> content = new HashMap<>();
    	
    	BufferedReader br = Files.newBufferedReader(Paths.get(gene_set_sum_file));
    	String header = br.readLine() + "\tp_Wilcoxon";
    	String line;
    	String set;
    	while((line=br.readLine())!=null) {
    		set = line.split("\t")[0];
    		content.put(set,line);
    	}
    	br.close();
    	
    	ValueComparator vc = new ValueComparator(set2pval);
    	TreeMap<String,Double> sorted_pvals = new TreeMap<>(vc);
    	sorted_pvals.putAll(set2pval);
    	BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
    	bw.write(header);
    	bw.newLine();
    	
    	for(String s: sorted_pvals.keySet()) {
    		bw.write(content.get(s)+"\t"+set2pval.get(s));
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
		return "Path2ExtendedSetSummary";
	}

	@Override
	protected String getOutfile() {
		return this.outfile;
	}

	@Override
	protected void saveExtraSettingsTo(NodeSettingsWO settings) {	
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings) throws InvalidSettingsException {
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings) throws InvalidSettingsException {
	}
}
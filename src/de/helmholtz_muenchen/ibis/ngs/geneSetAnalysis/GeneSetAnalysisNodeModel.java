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
import java.util.List;
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

import de.helmholtz_muenchen.ibis.utils.ngs.ContingencyTable;
import de.helmholtz_muenchen.ibis.ngs.lofsummary.Identifier;
import de.helmholtz_muenchen.ibis.ngs.lofsummary.MatrixSummary;
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
	static final String CFGKEY_GENE_SET_INFILE = "gene_set_infile";
		
	final SettingsModelString m_genesetin = new SettingsModelString(GeneSetAnalysisNodeModel.CFGKEY_GENE_SET_INFILE,"");

	protected static final NodeLogger logger = NodeLogger.getLogger(GeneSetAnalysisNodeModel.class);
	
	private String outfile;
	
    protected GeneSetAnalysisNodeModel() {
    	super();
    }

	@Override
	protected void performAnalysis(BufferedDataTable[] inData, ExecutionContext exec,
			HashMap<String, Double> gene2frequency, int pop_size, MatrixSummary ms)
					throws IOException, InvalidSettingsException {
		//get input
    	String matrix_file = inData[0].iterator().next().getCell(0).toString();
    	outfile = IO.replaceFileExtension(matrix_file, ".gene_set_analysis.tsv");
    	
    	// check/read gene set file
    	logger.debug("Read gene set file...");
    	String geneset_file = m_genesetin.getStringValue();
    	if(geneset_file.equals("") || Files.notExists(Paths.get(geneset_file))) {
    		throw new InvalidSettingsException("No gene set file specified!");
    	}

    	HashMap<String, ContingencyTable> gene2table = ms.toTables(new GeneSymbolIdentifier());
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
    	
    	HashMap<String, ContingencyTable> set2table = ms.toTables(new GeneSetIdentifier(genes2sets));
    	
    	//store input and results
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
    	
    	writeResults(outfile, set2genes,set2table, set2pval); 
	}

	class GeneSymbolIdentifier implements Identifier {

		@Override
		public List<String> getMappings(String identifier) {
			List<String> res = new ArrayList<String>();
			res.add(identifier.split("_")[2].toUpperCase());
			return res;
		}
	}
	
	class GeneSetIdentifier implements Identifier {

		HashMap <String, HashSet<String>> gene2sets;
		
		public GeneSetIdentifier(HashMap <String, HashSet<String>> gene2sets) {
			this.gene2sets = gene2sets;
		}
		
		@Override
		public List<String> getMappings(String identifier) {
			String gene = identifier.split("_")[2].toUpperCase();
			List<String> res = new ArrayList<String>();
			if(gene2sets.containsKey(gene)) {
				res.addAll(gene2sets.get(gene));
			}
			return res;
		}
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
	
//	private HashMap<String, HashSet<String>> readGeneSetSummary(String gene_set_sum_file) throws IOException {
//		HashMap<String, HashSet<String>> res = new HashMap<>();
//    	
//    	BufferedReader br = Files.newBufferedReader(Paths.get(gene_set_sum_file));
//		String line;
//		String [] fields;
//		HashSet<String> tmp;
//		while((line=br.readLine())!=null) {
//			fields = line.split("\t");
//			tmp = new HashSet<>();
//			for(String s: fields[SET_INDEX].split(",")) {
//				tmp.add(s.trim());
//			}
//			res.put(fields[SET_ID_INDEX],tmp);
//		}
//		return res;
//	}
	
    private void writeResults(String outfile, HashMap<String, HashSet<String>> set2genes, HashMap<String, ContingencyTable> set2table, HashMap<String, Double> set2pval) throws IOException {
//    	HashMap<String,String> content = new HashMap<>();
    	
//    	BufferedReader br = Files.newBufferedReader(Paths.get(gene_set_sum_file));
//    	String header = br.readLine() + "\tp_Wilcoxon";
    	String header = "set\tsize\tp_Wilcoxon";
//    	String line;
//    	String set;
//    	while((line=br.readLine())!=null) {
//    		set = line.split("\t")[0];
//    		content.put(set,line);
//    	}
//    	br.close();
    	
    	ValueComparator vc = new ValueComparator(set2pval);
    	TreeMap<String,Double> sorted_pvals = new TreeMap<>(vc);
    	sorted_pvals.putAll(set2pval);
    	BufferedWriter bw = Files.newBufferedWriter(Paths.get(outfile));
    	bw.write(header);
    	bw.newLine();
    	
    	for(String s: sorted_pvals.keySet()) {
    		bw.write(s+"\t"+set2genes.get(s).size()+"\t"+set2table.get(s).verticalToString()+"\t"+set2pval.get(s));
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
		m_genesetin.saveSettingsTo(settings);
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings) throws InvalidSettingsException {
		m_genesetin.loadSettingsFrom(settings);
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings) throws InvalidSettingsException {
		m_genesetin.validateSettings(settings);
	}
}
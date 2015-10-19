package de.helmholtz_muenchen.ibis.ngs.geneticBackgroundModel;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
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

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFFile;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFVariant;
import de.helmholtz_muenchen.ibis.utils.ngs.VEPAnnotationParser;

/**
 * This is the model implementation of GeneticBackgroundModel.
 * 
 *
 * @author Tim Jeske
 */
public class GeneticBackgroundModelNodeModel extends NodeModel {
	
	//configuration keys
	static final String CFGKEY_GTF_AFF = "gtf_aff";
	static final String CFGKEY_AC = "ac";
	static final String CFGKEY_AN = "an";
	
	//settings models
	static final String [] BASIS = new String[]{"genotypes","allele frequencies"};
	private final SettingsModelString m_gtf_aff = new SettingsModelString(CFGKEY_GTF_AFF, BASIS[0]);
	private final SettingsModelString m_ac = new SettingsModelString(CFGKEY_AC,"AC");
	private final SettingsModelString m_an = new SettingsModelString(CFGKEY_AN,"AN");
    
	private int vcf_index;
	
	//output col names
	public static final String OUT_COL1 = "Path2GeneticBackgroundModel";
	
	NodeLogger LOGGER = NodeLogger.getLogger(GeneticBackgroundModelNodeModel.class);
	
    /**
     * Constructor for the node model.
     */
    protected GeneticBackgroundModelNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String vcf_infile, outfile;
    	
    	vcf_infile = inData[0].iterator().next().getCell(vcf_index).toString();
    	
    	if(Files.notExists(Paths.get(vcf_infile))) {
    		throw new InvalidSettingsException("Input VCF file does not exist!");
    	}
    	
    	VCFFile vcf_it = new VCFFile(vcf_infile);
    	String vep_header = vcf_it.getInfoHeader(VEPAnnotationParser.INFO_ID);
    	if(vep_header == null) {
    		throw new InvalidSettingsException("No VEP annotations found!");
    	}
    	
    	VEPAnnotationParser parser = new VEPAnnotationParser(vep_header);
    	
    	HashMap<String, Double> gene_frequency;
    	
    	if(m_gtf_aff.getStringValue().equals(BASIS[0])) {
    		gene_frequency = fillGenotypes(vcf_it, parser);
    		outfile = IO.replaceFileExtension(vcf_infile, ".gene_model_gtf.tsv");
    	} else {
    		gene_frequency = fillAF(vcf_it, m_ac.getStringValue(), m_an.getStringValue(), parser);
    		outfile = IO.replaceFileExtension(vcf_infile, ".gene_model_aff.tsv");
    	}
    	
    	writeModel(outfile, gene_frequency);
    	
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

	private HashMap<String, Double> fillAF(VCFFile vcf_it, String ac_id, String an_id, VEPAnnotationParser parser) {
		HashMap<String, Double> result = new HashMap<>();
		VCFVariant var;
		String ac, an, csq;
		while(vcf_it.hasNext()) {
			var = vcf_it.next();
			ac = var.getInfoField(ac_id);
			an = var.getInfoField(an_id);
			csq = var.getInfoField(VEPAnnotationParser.INFO_ID);
			
			if(ac == null || an == null) {
				LOGGER.error("The INFO fields "+ac_id+ " and "+an_id+ " have not been found for variant on chr "+var.getChrom()+" at position "+var.getPos()+"!");
				continue;
			}
			
			if(csq == null) {
				LOGGER.error("No annotations have been found for variant on chr "+var.getChrom()+" at position "+var.getPos()+"!");
				continue;
			}
			
			HashMap<String, HashSet<Integer>> gene2allele_num = parser.getGene2AlleleNum(csq);
			String [] acs = ac.split(",");
			
			//compute frequency of being unaffected for each gene 
			double ac_abs, af, not_aff;
			
			for(String g: gene2allele_num.keySet()) {
				ac_abs = 0;
				for(int i: gene2allele_num.get(g)) {
					ac_abs += Integer.parseInt(acs[i-1]);
				}
				
				af = (double)ac_abs/Double.parseDouble(an);
				not_aff = Math.pow(1.0-af,2);
				if(result.containsKey(g)) {
					not_aff = result.get(g) * not_aff;
				}
				result.put(g, not_aff);
			}
			
		}
		
		for(String g: result.keySet()) {
			result.put(g, 1.0 - result.get(g));
		}
		
		return result;
	}

	private HashMap<String, Double> fillGenotypes(VCFFile vcf_it, VEPAnnotationParser parser) {
		HashMap<String, Double> result = new HashMap<>();
		HashMap<String, HashSet<String>> unaffected_samples, affected_samples;
		unaffected_samples = new HashMap<>();
		affected_samples = new HashMap<>();
		VCFVariant var;
		String csq;
		while(vcf_it.hasNext()) {
			var = vcf_it.next();
			csq = var.getInfoField(VEPAnnotationParser.INFO_ID);
			
			if(csq == null) {
				LOGGER.error("No annotations have been found for variant on chr "+var.getChrom()+" at position "+var.getPos()+"!");
				continue;
			}
			
			HashMap<String, HashSet<Integer>> gene2allele_num = parser.getGene2AlleleNum(csq);
			
			HashSet<Integer> alt_alleles;
			HashSet<String> aff_samp, unaff_samp;
			
			for(String gene: gene2allele_num.keySet()) {
				alt_alleles = gene2allele_num.get(gene);
				aff_samp = var.getAffectedSamples(alt_alleles);
				unaff_samp = var.getUnaffectedSamples(alt_alleles);
				if(unaffected_samples.containsKey(gene)) {
					//unaffected samples contains all samples that carry only non-LOF alleles in a gene
					unaffected_samples.get(gene).retainAll(unaff_samp);
				} else {
					unaffected_samples.put(gene,unaff_samp);
				}
				
				if(affected_samples.containsKey(gene)) {
					affected_samples.get(gene).addAll(aff_samp);
				} else {
					affected_samples.put(gene, aff_samp);
				}
			}
		}
		
		HashSet<String> called_samples;
		for(String gene: affected_samples.keySet()) {
			called_samples = new HashSet<>();
			called_samples.addAll(unaffected_samples.get(gene));
			called_samples.addAll(affected_samples.get(gene));
			result.put(gene, (double)affected_samples.get(gene).size()/(double)called_samples.size());
		}
	
		return result;
	}
	
    private void writeModel(String outfile, HashMap<String, Double> gene_frequency) {
		
    	Charset c = Charset.forName("UTF-8");
    	BufferedWriter bw;
    	try {
			bw = Files.newBufferedWriter(Paths.get(outfile), c);
			for(String s: gene_frequency.keySet()) {
				bw.write(s+"\t"+gene_frequency.get(s));
				bw.newLine();
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

    	vcf_index = -1;
    	
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
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_gtf_aff.saveSettingsTo(settings);
    	m_ac.saveSettingsTo(settings);
    	m_an.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_gtf_aff.loadSettingsFrom(settings);
    	m_ac.loadSettingsFrom(settings);
    	m_an.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_gtf_aff.validateSettings(settings);
    	m_ac.validateSettings(settings);
    	m_an.validateSettings(settings);
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

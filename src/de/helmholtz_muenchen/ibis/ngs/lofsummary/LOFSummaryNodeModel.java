package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

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
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.AnnotationParser;
import de.helmholtz_muenchen.ibis.utils.ngs.VEPAnnotationParser;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;
import de.helmholtz_muenchen.ibis.utils.ngs.VCFFile;

/**
 * This is the model implementation of LOFSummarizer.
 * 
 *
 * @author Tim Jeske
 */
public class LOFSummaryNodeModel extends NodeModel {
    
	// the logger instance
    protected static final NodeLogger logger = NodeLogger.getLogger(LOFSummaryNodeModel.class);
	
	static final String CFGKEY_CDS_INFILE = "cds_infile";
	static final String CFGKEY_PED_INFILE = "ped_infile";
	static final String CFGKEY_GENE_SET_INFILE = "gene_set_infile";
	
	final SettingsModelString m_cdsin = new SettingsModelString(LOFSummaryNodeModel.CFGKEY_CDS_INFILE,"");
	final SettingsModelString m_pedin = new SettingsModelString(LOFSummaryNodeModel.CFGKEY_PED_INFILE,"");
	final SettingsModelString m_genesetin = new SettingsModelString(LOFSummaryNodeModel.CFGKEY_GENE_SET_INFILE,"");
	
	//selected annotation
    static final String CFGKEY_ANNOTATION="annotation";
    static final String[] ANNOTATIONS_AVAILABLE={"VEP"};
    final SettingsModelString m_annotation = new SettingsModelString(CFGKEY_ANNOTATION, "");
	
	//output
	public static final String OUT_COL1 = "Path2Gene_Summary";
	
	private int vcf_index;
	
    /**
     * Constructor for the node model.
     */
    protected LOFSummaryNodeModel() {
    
    	super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	String vcf_infile = inData[0].iterator().next().getCell(vcf_index).toString();
    	
    	if(Files.notExists(Paths.get(vcf_infile))) {
    		throw new InvalidSettingsException("Input VCF file does not exist!");
    	}
    	
    	String infile_warning = CheckUtils.checkSourceFile(vcf_infile);
    	if(infile_warning != null) {
    		setWarningMessage(infile_warning);
    	}

    	String cds_file = m_cdsin.getStringValue();
    	if(cds_file.equals("") || Files.notExists(Paths.get(cds_file))) {
    		cds_file = null;
    		throw new InvalidSettingsException("No CDS or GTF file specified!");
    	}
    	
    	String ped_file = m_pedin.getStringValue();
    	
    	if(ped_file.equals("") || Files.notExists(Paths.get(ped_file))) {
    		setWarningMessage("No PED file specified!");
    		ped_file = null;
    	}
    	
    	String geneset_file = m_genesetin.getStringValue();
    	if(geneset_file.equals("") || Files.notExists(Paths.get(geneset_file))) {
    		setWarningMessage("No gene set file specified!");
    		geneset_file = null;
    	}

    	String var_outfile = IO.replaceFileExtension(vcf_infile, ".variant_summary.tsv");
    	String gene_outfile = IO.replaceFileExtension(vcf_infile, ".gene_summary.tsv");
    	String sample_outfile = IO.replaceFileExtension(vcf_infile, ".sample_summary.tsv");
    	String enrichment_outfile = IO.replaceFileExtension(vcf_infile, ".gene_set_summary.tsv");
    	
    	VCFFile vcf = new VCFFile(vcf_infile);
    	AnnotationParser parser = null;
    	String annotation = m_annotation.getStringValue();
    	if(annotation.equals("VEP")) {
    		parser = new VEPAnnotationParser(vcf.getInfoHeader(VEPAnnotationParser.ANN_ID));
    	}
    	
    	logger.debug("Read CDS/GTF file...");
    	HashMap<String, Gene> geneid2gene = new HashMap<>();
    	if(cds_file.endsWith("fa") || cds_file.endsWith("fasta")) {
    		geneid2gene = this.readCDSFile(cds_file);
    	} else if(cds_file.endsWith("gtf")) {
    		geneid2gene = this.readGTFFile(cds_file);
    	}
    	
    	logger.debug("Read PED file...");
    	HashMap<String, Boolean> sampleid2affected = readPEDFile(ped_file);
    	
    	logger.debug("Read gene set file...");
    	HashMap<String, HashSet<String>> set2genes = readGeneSetFile(geneset_file);
    	
    	//parallel execution
    	final AnnotationParser myparser = parser;
    	final HashMap<String, Gene> mygene2gene = geneid2gene;
    	
    	Thread t1 = new Thread() {
    		public void run() {
    			try {
    				logger.info("Generate variant summary...");
					LOFSummarizer.getVarSum(vcf, myparser, mygene2gene, var_outfile);
					logger.info("Variant summary ready!");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
    		}
    	};
    	
    	Thread t2 = new Thread() {
    		public void run() {
    	    	try {
    	    		logger.info("Generate gene summary...");
					LOFSummarizer.getGeneSum(new VCFFile(vcf_infile), myparser, mygene2gene, sampleid2affected, gene_outfile);
					logger.info("Gene summary ready!");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
    		}
    	};
    	
    	Thread t3 = new Thread() {
    		public void run() {
    	    	try {
    	    		logger.info("Generate sample summary...");
    	    		LOFSummarizer.getSampleSum(new VCFFile(vcf_infile), myparser, sample_outfile);
					logger.info("Sample summary ready!");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
    		}
    	};
    	Thread t4 = new Thread() {
    		public void run() {
    	    	try {
    	    		logger.info("Generate gene set summary...");
    	    		LOFSummarizer.getEnrichmentSum(new VCFFile(vcf_infile), myparser, set2genes, sampleid2affected, enrichment_outfile);
					logger.info("Gene set summary ready!");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
    		}
    	};
    	
    	t1.start();
    	t2.start();
    	t3.start();
    	t4.start();
    	
    	t1.join();
    	t2.join();
    	t3.join();
    	t4.join();
    	
    	//Create Output Table
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(gene_outfile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
    	
        return new BufferedDataTable[]{outTable};
    }
    
    /**
	 * reads CDS file and returns map with gene ids as keys and Gene objects as values
	 * @throws IOException
	 */
	private HashMap<String, Gene> readCDSFile(String cds_file) throws IOException {
		HashMap<String, Gene> geneid2gene = new HashMap<>();
		String [] fields;
		String transcript_id, gene_id;
		
		FileInputStream inputStream = null;
		Scanner sc = null;
		try {
		    inputStream = new FileInputStream(cds_file);
		    sc = new Scanner(inputStream, "UTF-8");
		    while (sc.hasNextLine()) {
		        String line = sc.nextLine();
		        if(line.startsWith(">") && line.contains("transcript_biotype:protein_coding")) {
		        	fields = line.split("\\s");
		        	transcript_id = fields[0].replaceFirst(">","");
		        	gene_id = fields[3].split(":")[1];
		        	
		        	if(geneid2gene.containsKey(gene_id)) {
		        		geneid2gene.get(gene_id).addTranscript(transcript_id);
		        	} else {
		        		Gene g = new Gene(gene_id);
		        		g.addTranscript(transcript_id);
		        		geneid2gene.put(gene_id,g);
		        	}
		        }
		    }
		    // note that Scanner suppresses exceptions
		    if (sc.ioException() != null) {
		        throw sc.ioException();
		    }
		} finally {
		    if (inputStream != null) {
		        inputStream.close();
		    }
		    if (sc != null) {
		        sc.close();
		    }
		}
		return geneid2gene;
	}
	
    /**
	 * reads GTF file and returns map with gene ids as keys and Gene objects as values
	 * @throws IOException
	 */
	private HashMap<String, Gene> readGTFFile(String gtf_file) throws IOException {
		HashMap<String, Gene> geneid2gene = new HashMap<>();
		String [] fields, info_fields;
		String transcript_id = "";
		String gene_id = "";
		String gene_symbol = "";
		
		FileInputStream inputStream = null;
		Scanner sc = null;
		try {
		    inputStream = new FileInputStream(gtf_file);
		    sc = new Scanner(inputStream, "UTF-8");
		    while (sc.hasNextLine()) {
		        String line = sc.nextLine();
		        if(line.startsWith("#")) continue;
		        fields = line.split("\t");
		        if(fields[1].equals("protein_coding")) {
		        	if(fields[8].contains("gene_id") && fields[8].contains("transcript_id") && fields[8].contains("gene_name")) {
						info_fields = fields[8].split(";");
						for (String i : info_fields) {
							i = i.trim();
							if (i.startsWith("gene_id")) {
								gene_id = i.split("\\s")[1].replaceAll("\"", "");
							} else if (i.startsWith("transcript_id")) {
								transcript_id = i.split("\\s")[1].replaceAll("\"", "");
							} else if (i.startsWith("gene_name")) {
								gene_symbol = i.split("\\s")[1].replaceAll("\"", "");
							}
						}

						if (geneid2gene.containsKey(gene_id)) {
							geneid2gene.get(gene_id).addTranscript(transcript_id);
						} else {
							Gene g = new Gene(gene_id, gene_symbol);
							g.addTranscript(transcript_id);
							geneid2gene.put(gene_id, g);
						}
		        	}
		        }
		    }
		    // note that Scanner suppresses exceptions
		    if (sc.ioException() != null) {
		        throw sc.ioException();
		    }
		} finally {
		    if (inputStream != null) {
		        inputStream.close();
		    }
		    if (sc != null) {
		        sc.close();
		    }
		}
		return geneid2gene;
	}
	
	private HashMap<String, Boolean> readPEDFile(String ped_file) throws IOException {
		HashMap<String, Boolean> res = new HashMap<>();
		
		BufferedReader br = Files.newBufferedReader(Paths.get(ped_file));
		String line;
		String[] fields;
		while ((line = br.readLine()) != null) {
			if (line.startsWith("Family ID"))
				continue;

			fields = line.split("\t");
			if (fields[5].equals("2")) {
				res.put(fields[1], true);
			} else {
				res.put(fields[1], false);
			}
		}
		br.close();
		
		return res;
	}
	
	private HashMap<String, HashSet<String>> readGeneSetFile(String file) throws IOException {
		HashMap<String, HashSet<String>> set2genes = new HashMap<>();
		
		BufferedReader br = Files.newBufferedReader(Paths.get(file));
		String line;
		String [] fields;
		
		while((line=br.readLine())!=null) {
			line = line.trim();
			fields = line.split("\t");
			HashSet<String> tmp = new HashSet<>();
			for(int i = 2; i < fields.length; i++) {
				tmp.add(fields[i]);
			}
			set2genes.put(fields[0], tmp);
		}
		
		br.close();
		
		return set2genes;
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
    	m_cdsin.saveSettingsTo(settings);
    	m_pedin.saveSettingsTo(settings);
    	m_annotation.saveSettingsTo(settings);
    	m_genesetin.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_cdsin.loadSettingsFrom(settings);
    	m_pedin.loadSettingsFrom(settings);
    	m_annotation.loadSettingsFrom(settings);
    	m_genesetin.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_cdsin.validateSettings(settings);
    	m_pedin.validateSettings(settings);
    	m_annotation.validateSettings(settings);
    	m_genesetin.validateSettings(settings);
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
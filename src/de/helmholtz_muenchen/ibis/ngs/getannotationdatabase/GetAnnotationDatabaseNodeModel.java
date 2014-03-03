package de.helmholtz_muenchen.ibis.ngs.getannotationdatabase;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
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

import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of GetAnnotationDatabase.
 * 
 * @author Sebastian Kopetzky
 */
public class GetAnnotationDatabaseNodeModel extends NodeModel {
	
	public static final String CFGKEY_INSTALLPATH = "installpath";
	public static final String CFGKEY_ORGANISM = "organism";
	public static final String CFGKEY_DATABASE = "database";
	public static final String CFGKEY_ANNOTATIONTYPE = "annoType";
	public static final String CFGKEY_OUTNAME = "outname";

	private final SettingsModelString m_installpath = new SettingsModelString(CFGKEY_INSTALLPATH,"");
	private final SettingsModelString m_organism = new SettingsModelString(CFGKEY_ORGANISM,"");
	private final SettingsModelString m_database = new SettingsModelString(CFGKEY_DATABASE,"");
	private final SettingsModelString m_annotationtype = new SettingsModelString(CFGKEY_ANNOTATIONTYPE,"");
	private final SettingsModelString m_outname = new SettingsModelString(CFGKEY_OUTNAME,"");
	private static final NodeLogger LOGGER = NodeLogger.getLogger(GetAnnotationDatabaseNodeModel.class);
	
    /**
     * Constructor for the node model.
     */
    protected GetAnnotationDatabaseNodeModel() {
    	
        super(0, 0);
        
        m_outname.setStringValue("humandb");
        m_database.setStringValue("hg18");
        m_outname.setEnabled(false);
        m_database.setEnabled(false);
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	/**Initialize logfile**/
    	String logfile = m_installpath.getStringValue()+"/logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(ShowOutput.getNodeStartTime("GetAnnotationDatabase"));
    	ShowOutput.writeLogFile(logBuffer);
    	/**end initializing logfile**/
    	
    	ArrayList<String> command = new ArrayList<String>();
    	String installpath = m_installpath.getStringValue() + "/";
    	String program = "annotate_variation.pl ";
    	String database = m_database.getStringValue();
    	String annoType = m_annotationtype.getStringValue();
    	String outname = m_outname.getStringValue();
    	//String organism = m_organism.getStringValue();
    	
    	if(annoType.equals("[UCSC:gene-based] RefSeq transcript annotations and mRNA sequences in FASTA formats")) {
    		annoType = "gene";
    	} else if(annoType.equals("[UCSC:gene-based] UCSC Gene annotation (more comprehensive than RefSeq annotation) and mRNA sequences in FASTA format")) {
    		annoType = "knownGene";
    	} else if(annoType.equals("[UCSC:gene-based] Ensembl Gene annotation (more comprehensive than RefSeq annotation) and mRNA sequences in FASTA format")) {
    		annoType = "ensGene";
    	} else if(annoType.equals("[UCSC:gene-based] GENCODE manual annotation for hg18 (human)")) {
    		annoType = "wgEncodeGencodeManualV3";
    	} else if(annoType.equals("[UCSC:gene-based] GENCODE manual annotation for hg19 (human)")) {
    		annoType = "wgEncodeGencodeManualV4";
    	} else if(annoType.equals("[UCSC:region-based] Approximate location of bands seen on Giemsa-stained chromosomes")) {
    		annoType = "band";
    	} else if(annoType.equals("[UCSC:region-based] Transcription factor binding sites conserved in the human/mouse/rat alignment, based on transfac Matrix Database")) {
    		annoType = "tfbs";
    	} else if(annoType.equals("[UCSC:region-based] snoRNA and miRNA annotations")) {
    		annoType = "mirna";
    	} else if(annoType.equals("[UCSC:region-based] TargetScan generated miRNA target site predictions")) {
    		annoType = "mirnatarget";
    	} else if(annoType.equals("[UCSC:region-based] Segmental duplications in genome")) {
    		annoType = "segdup";
    	} else if(annoType.equals("[UCSC:region-based] Conserved elements produced by the phastCons program based on a whole-genome alignment of vertebrates")) {
    		annoType = "mce28way";
    	} else if(annoType.equals("[UCSC:region-based] Conserved functional RNA, through RNA secondary structure predictions made with the EvoFold program")) {
    		annoType = "evofold";
    	} else if(annoType.equals("[UCSC:region-based] Database of Genomic Variants, which contains annotations for reported structural variations")) {
    		annoType = "dgv";
    	} else if(annoType.equals("[UCSC:region-based] Canonical UCSC genes that have been associated with identifiers in the OMIM database")) {
    		annoType = "omimgene";
    	} else if(annoType.equals("[UCSC:region-based] Published GWAS results on diverse human diseases")) {
    		annoType = "gwascatalog";
    	} else if(annoType.equals("[1000 Genomes Project] Pilot 1 allele frequency data on the CEU, YRI and JPTCHB populations, updated on April 2009")) {
    		annoType = "1000g";
    	} else if(annoType.equals("[1000 Genomes Project] Pilot 1 allele frequency data on the CEU, YRI and JPTCHB populations, updated on March 2010")) {
    		annoType = "1000g2010";
    	} else if(annoType.equals("[1000 Genomes Project] Pilot 1 allele frequency data on the CEU, YRI and JPTCHB populations, updated on July 2010")) {
    		annoType = "1000g2010jul";
    	} else if(annoType.equals("[1000 Genomes Project] Full phase 1 project with 629 subjects form diverse populations using August 2010, updated on November 2010")) {
    		annoType = "1000g2010nov";
    	} else if(annoType.equals("[1000 Genomes Project] PHASE 1 low-coverage data on 1094 subjects using November 2010 alignments, updated on May 2011")) {
    		annoType = "1000g2011may";
    	} else if(annoType.equals("[1000 Genomes Project] variant calls from 1092 samples for SNPs, short indels, updated on February 2012")) {
    		annoType = "1000g2012feb";
    	} else if(annoType.equals("[dbSNP] Version 132")) {
    		annoType = "snp132";
    	} else if(annoType.equals("[dbSNP] Version 131")) {
    		annoType = "snp131";
    	} else if(annoType.equals("[dbSNP] Version 130")) {
    		annoType = "snp130";
    	} else if(annoType.equals("[dbSNP] Version 129")) {
    		annoType = "snp129";
    	} else if(annoType.equals("[dbSNP] version 125")) {
    		annoType = "snp125";
    	} else if(annoType.equals("[SIFT] Chr, start, end, reference allele, observed allele, SIFT score, reference amino acid, observed amino acids, updated on March 2011")) {
    		annoType = "avsift";
    	} else if(annoType.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (PolyPhen2)")) {
    		annoType = "ljb_pp2 -webfrom annovar";
    	} else if(annoType.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (using LJBSIFT scores)")) {
    		annoType = "ljb_sift -webfrom annovar";
    	} else if(annoType.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (MutationTaster)")) {
    		annoType = "ljb_mt -webfrom annovar";
    	} else if(annoType.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (PhyloP conservation score)")) {
    		annoType = "ljb_phylop -webfrom annovar";
    	} else if(annoType.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (LRT)")) {
    		annoType = "ljb_lrt -webfrom annovar";
    	} else if(annoType.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (GERP++ score for exonic variants)")) {
    		annoType = "ljb_gerp++ -webfrom annovar";
    	} else if(annoType.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (all scores above from LJB database)")) {
    		annoType = "ljb_all -webfrom annovar";
    	} else if(annoType.equals("[ESP] 5400 NHLBI exomes (European Americans)")) {
    		annoType = "esp5400_ea -webfrom annovar";
    	} else if(annoType.equals("[ESP] 5400 NHLBI exomes (African Americans)")) {
    		annoType = "esp5400_aa -webfrom annovar";
    	} else if(annoType.equals("[ESP] 5400 NHLBI exomes (all ethnicity)")) {
    		annoType = "esp5400_all -webfrom annovar";
    	}


    	command.add(installpath+program);
    	command.add("-buildver "+database);
    	command.add("-downdb " + annoType);
    	command.add(installpath+outname);  	
    	
    	/**
    	 * Execute
    	 */
    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);
    	logBuffer.append(ShowOutput.getNodeEndTime());
    	ShowOutput.writeLogFile(logBuffer);
    	
        return new BufferedDataTable[]{};
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
    	
        return new DataTableSpec[]{};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_annotationtype.saveSettingsTo(settings);
    	m_database.saveSettingsTo(settings);
    	m_installpath.saveSettingsTo(settings);
    	m_organism.saveSettingsTo(settings);
    	m_outname.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_annotationtype.loadSettingsFrom(settings);
    	m_database.loadSettingsFrom(settings);
    	m_installpath.loadSettingsFrom(settings);
    	m_organism.loadSettingsFrom(settings);
    	m_outname.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_annotationtype.loadSettingsFrom(settings);
    	m_database.loadSettingsFrom(settings);
    	m_installpath.loadSettingsFrom(settings);
    	m_organism.loadSettingsFrom(settings);
    	m_outname.loadSettingsFrom(settings);
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


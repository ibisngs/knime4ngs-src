package de.helmholtz_muenchen.ibis.ngs.kggseq;

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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of KGGSeq.
 * 
 *
 * @author 
 */
public class KGGSeqNodeModel extends NodeModel {
    
	
	/**
	 * Logger
	 */
	private static final NodeLogger LOGGER = NodeLogger.getLogger(KGGSeqNodeModel.class);
	
	
	
	/**
	 * Config Keys
	 */
	
//	Steps according to http://statgenpro.psychiatry.hku.hk/limx/kggseq/doc/Denovo.htm
	
//	Step 1
	public static final String CFGKEY_KGGSEQ_PATH 			= "KGGSEQ_PATH";
	public static final String CFGKEY_BUILDVER			  	= "BUILDVER";
	public static final String CFGKEY_INFILE			  	= "INFILE";
	public static final String CFGKEY_PEDFILE 			   	= "PEDFILE";
	public static final String CFGKEY_COMPOSITE_SUBJECT_ID 	= "COMPOSITE_SUBJECT_ID";
	public static final String CFGKEY_OUTPREFIX				= "OUTPREFIX";
	public static final String CFGKEY_OUTFORMAT				= "OUTFORMAT";
	public static final String CFGKEY_SEQ_QUAL				= "SEQ_QUAL";
	public static final String CFGKEY_SEQ_MQ				= "SEQ_MQ";
	public static final String CFGKEY_SEQ_SB				= "SEQ_SB";
	public static final String CFGKEY_GTY_QUAL				= "GTY_QUAL";
	public static final String CFGKEY_GTY_DP				= "GTY_DP";
	public static final String CFGKEY_GTY_SEC_PL			= "GTY_SEC_PL";
	public static final String CFGKEY_GTY_AF_REF			= "GTY_AF_REF";
	public static final String CFGKEY_GTY_AF_HET			= "GTY_AF_HET";
	public static final String CFGKEY_GTY_AF_ALT			= "GTY_AF_ALT";

//	Step 2
	public static final String CFGKEY_GENOTYPE_FILTER		= "GENOTYPE_FILTER";
	public static final String CFGKEY_IGNORE_HOMO			= "IGNORE_HOMO";
	
//  Step 3
	public static final String CFGKEY_GENE_FEATURES			= "GENE_FEATURES";
	
//	Step 4
	public static final String CFGKEY_FILTER_COMMON			= "FILTER_COMMON";
	
//	Step 5
	public static final String CFGKEY_DISEASE_CAUSING_PRED	= "DISEASE_CAUSING_PRED";

//	Step 6
	public static final String CFGKEY_OMIM_ANNO				= "OMIM_ANNO";
	
//	Step 7
	public static final String CFGKEY_CANDIDATE_PPI			= "CANDIDATE_PPI";
	public static final String CFGKEY_CANDIDATE_GENES		= "CANDIDATE_GENES";
	
//	Step 8
	public static final String CFGKEY_CANDIDATE_PATHWAYS	= "CANDIDATE_PATHWAYS";
	
//	Step 9
	public static final String CFGKEY_PUBMED	= "PUBMED";
	
	
	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_KGGSEQ = new SettingsModelString(KGGSeqNodeModel.CFGKEY_KGGSEQ_PATH, "");
    private final SettingsModelString m_BUILDVER = new SettingsModelString(KGGSeqNodeModel.CFGKEY_BUILDVER, "hg19");
    private final SettingsModelString m_INFILE = new SettingsModelString(KGGSeqNodeModel.CFGKEY_INFILE, "");
    private final SettingsModelString m_PEDFILE = new SettingsModelString(KGGSeqNodeModel.CFGKEY_PEDFILE, "");
    private final SettingsModelBoolean m_COMPOSITESUBJECTID = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_COMPOSITE_SUBJECT_ID, false);
    private final SettingsModelString m_OUTPREFIX = new SettingsModelString(KGGSeqNodeModel.CFGKEY_OUTPREFIX, "");
    private final SettingsModelString m_OUTFORMAT = new SettingsModelString(KGGSeqNodeModel.CFGKEY_OUTFORMAT, "excel");
   
    private final SettingsModelDoubleBounded m_SEQ_QUAL = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_SEQ_QUAL,50,0,Double.MAX_VALUE);
    private final SettingsModelDoubleBounded m_SEQ_MQ = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_SEQ_MQ,20,0,Double.MAX_VALUE);
    private final SettingsModelDoubleBounded m_SEQ_SB = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_SEQ_SB,-10,-Double.MAX_VALUE,Double.MAX_VALUE);
    private final SettingsModelDoubleBounded m_GTY_QUAL = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_QUAL,20,0,Double.MAX_VALUE);
    private final SettingsModelDoubleBounded m_GTY_DP = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_DP,8,0,Double.MAX_VALUE);
    private final SettingsModelIntegerBounded m_GTY_SEC_PL = new SettingsModelIntegerBounded(KGGSeqNodeModel.CFGKEY_GTY_SEC_PL,20,0,Integer.MAX_VALUE);
    private final SettingsModelDoubleBounded m_GTY_AF_REF = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_AF_REF,0.05,0,Double.MAX_VALUE);
    private final SettingsModelDoubleBounded m_GTY_AF_HET = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_AF_HET,0.25,0,Double.MAX_VALUE);
    private final SettingsModelDoubleBounded m_GTY_AF_ALT = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_AF_ALT,0.75,0,Double.MAX_VALUE);

    private final SettingsModelOptionalString m_GENOTYPE_FILTER = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_GENOTYPE_FILTER, "4",true);
    private final SettingsModelBoolean m_IGNORE_HOMO = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_IGNORE_HOMO, true);
    private final SettingsModelOptionalString m_GENE_FEATURES = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_GENE_FEATURES, "0,1,2,3,4,5",true);
    private final SettingsModelBoolean m_FILTER_COMMON = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_FILTER_COMMON, true);   
    private final SettingsModelBoolean m_DISEASE_CAUSING_PRED = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_DISEASE_CAUSING_PRED, true);
    private final SettingsModelBoolean m_OMIM_ANNO = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_OMIM_ANNO, true);
    private final SettingsModelBoolean m_CANDIDATE_PPI = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_CANDIDATE_PPI, true);
    private final SettingsModelString m_CANDIDATE_GENES = new SettingsModelString(KGGSeqNodeModel.CFGKEY_CANDIDATE_GENES, "");
    private final SettingsModelBoolean m_CANDIDATE_PATHWAYS = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_CANDIDATE_PATHWAYS, true);
    private final SettingsModelOptionalString m_PUBMED = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_PUBMED,"",false);

	
    /**
     * Constructor for the node model.
     */
    protected KGGSeqNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	StringBuffer stdErr = new StringBuffer();
    	StringBuffer stdOut = new StringBuffer();
    	
    	Executor.executeCommand(createCommand(),exec,LOGGER,stdOut,stdErr);
    	
        return inData;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO: generated method stub
    }

    
    /**
     * Creates the execution command of the KGGSeq Node
     * @return
     */
    private String[] createCommand(){
    	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	command.add("java -jar");
    	command.add(m_KGGSEQ.getStringValue());
    	
    	command.add("--buildver "+m_BUILDVER.getStringValue());
    	command.add("--vcf-file "+m_INFILE.getStringValue());
    	command.add("--ped-file "+m_PEDFILE.getStringValue());
    	
    	if(m_COMPOSITESUBJECTID.getBooleanValue()){
    		command.add("--composite-subject-id");
    	}
    	
    	command.add("--out "+m_OUTPREFIX.getStringValue());
    	command.add("--"+m_OUTFORMAT.getStringValue());
    	command.add("--seq-qual "+m_SEQ_QUAL.getDoubleValue());
    	command.add("--seq-mq "+m_SEQ_MQ.getDoubleValue());
    	command.add("--seq-sb "+m_SEQ_SB.getDoubleValue());
    	command.add("--gty-qual "+m_GTY_QUAL.getDoubleValue());
    	command.add("--gty-dp "+m_GTY_DP.getDoubleValue());
    	command.add("--gty-sec-pl "+m_GTY_SEC_PL.getIntValue());
    	command.add("--gty-af-ref "+m_GTY_AF_REF.getDoubleValue());
    	command.add("--gty-af-het "+m_GTY_AF_HET.getDoubleValue());
    	command.add("--gty-af-alt "+m_GTY_AF_ALT.getDoubleValue());
    
    	if(m_GENOTYPE_FILTER.isActive() && !m_GENOTYPE_FILTER.getStringValue().equals("")){
    		command.add("--genotype-filter 7,"+m_GENOTYPE_FILTER.getStringValue());
    	}
    	
    	if(m_IGNORE_HOMO.getBooleanValue()){
    		command.add("--ignore-homo");
    	}
    	
    	if(m_GENE_FEATURES.isActive()){
    		command.add("--db-gene refgene --gene-feature-in "+m_GENE_FEATURES.getStringValue());
    	}
    	
    	if(m_FILTER_COMMON.getBooleanValue()){
    		command.add("--db-filter hg19_1kg201204,hg19_dbsnp138,hg19_ESP6500AA,hg19_ESP6500EA  --rare-allele-freq 0.05");
    	}
    	
    	if(m_DISEASE_CAUSING_PRED.getBooleanValue()){
    		command.add("--db-score dbnsfp --mendel-causing-predict all --filter-nondisease-variant");
    	}
    	
    	if(m_OMIM_ANNO.getBooleanValue()){
    		command.add("--genome-annot --omim-annot");
    	}
    	
    	if(m_CANDIDATE_GENES.isEnabled()){
    		command.add("--candi-file "+m_CANDIDATE_GENES.getStringValue());
    	}
    	
    	if(m_CANDIDATE_PPI.getBooleanValue()){
    		command.add("--ppi-annot string --ppi-depth 1");
    	}
    	
    	if(m_CANDIDATE_PATHWAYS.getBooleanValue()){
    		command.add("--pathway-annot cura");
    	}
    	
    	if(m_PUBMED.isActive()){
    		command.add("--pubmed-mining "+m_PUBMED.getStringValue());
    	}
    	
    	command.add("--no-resource-check");
    	
    	System.out.println(StringUtils.join(command, " "));
    	
    	return new String[]{StringUtils.join(command, " ")};
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

        // TODO: generated method stub
        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_BUILDVER.saveSettingsTo(settings);
         m_CANDIDATE_GENES.saveSettingsTo(settings);
         m_CANDIDATE_PATHWAYS.saveSettingsTo(settings);
         m_CANDIDATE_PPI.saveSettingsTo(settings);
         m_COMPOSITESUBJECTID.saveSettingsTo(settings);
         m_DISEASE_CAUSING_PRED.saveSettingsTo(settings);
         m_FILTER_COMMON.saveSettingsTo(settings);
         m_GENE_FEATURES.saveSettingsTo(settings);
         m_GENOTYPE_FILTER.saveSettingsTo(settings);
         m_GTY_AF_ALT.saveSettingsTo(settings);
         m_GTY_AF_HET.saveSettingsTo(settings);
         m_GTY_AF_REF.saveSettingsTo(settings);
         m_GTY_DP.saveSettingsTo(settings);
         m_GTY_QUAL.saveSettingsTo(settings);
         m_GTY_SEC_PL.saveSettingsTo(settings);
         m_IGNORE_HOMO.saveSettingsTo(settings);
         m_INFILE.saveSettingsTo(settings);
         m_KGGSEQ.saveSettingsTo(settings);
         m_OMIM_ANNO.saveSettingsTo(settings);
         m_OUTFORMAT.saveSettingsTo(settings);
         m_OUTPREFIX.saveSettingsTo(settings);
         m_PEDFILE.saveSettingsTo(settings);
         m_PUBMED.saveSettingsTo(settings);
         m_SEQ_MQ.saveSettingsTo(settings);
         m_SEQ_QUAL.saveSettingsTo(settings);
         m_SEQ_SB.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_BUILDVER.loadSettingsFrom(settings);
        m_CANDIDATE_GENES.loadSettingsFrom(settings);
        m_CANDIDATE_PATHWAYS.loadSettingsFrom(settings);
        m_CANDIDATE_PPI.loadSettingsFrom(settings);
        m_COMPOSITESUBJECTID.loadSettingsFrom(settings);
        m_DISEASE_CAUSING_PRED.loadSettingsFrom(settings);
        m_FILTER_COMMON.loadSettingsFrom(settings);
        m_GENE_FEATURES.loadSettingsFrom(settings);
        m_GENOTYPE_FILTER.loadSettingsFrom(settings);
        m_GTY_AF_ALT.loadSettingsFrom(settings);
        m_GTY_AF_HET.loadSettingsFrom(settings);
        m_GTY_AF_REF.loadSettingsFrom(settings);
        m_GTY_DP.loadSettingsFrom(settings);
        m_GTY_QUAL.loadSettingsFrom(settings);
        m_GTY_SEC_PL.loadSettingsFrom(settings);
        m_IGNORE_HOMO.loadSettingsFrom(settings);
        m_INFILE.loadSettingsFrom(settings);
        m_KGGSEQ.loadSettingsFrom(settings);
        m_OMIM_ANNO.loadSettingsFrom(settings);
        m_OUTFORMAT.loadSettingsFrom(settings);
        m_OUTPREFIX.loadSettingsFrom(settings);
        m_PEDFILE.loadSettingsFrom(settings);
        m_PUBMED.loadSettingsFrom(settings);
        m_SEQ_MQ.loadSettingsFrom(settings);
        m_SEQ_QUAL.loadSettingsFrom(settings);
        m_SEQ_SB.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_BUILDVER.loadSettingsFrom(settings);
        m_CANDIDATE_GENES.loadSettingsFrom(settings);
        m_CANDIDATE_PATHWAYS.loadSettingsFrom(settings);
        m_CANDIDATE_PPI.loadSettingsFrom(settings);
        m_COMPOSITESUBJECTID.loadSettingsFrom(settings);
        m_DISEASE_CAUSING_PRED.loadSettingsFrom(settings);
        m_FILTER_COMMON.loadSettingsFrom(settings);
        m_GENE_FEATURES.loadSettingsFrom(settings);
        m_GENOTYPE_FILTER.loadSettingsFrom(settings);
        m_GTY_AF_ALT.loadSettingsFrom(settings);
        m_GTY_AF_HET.loadSettingsFrom(settings);
        m_GTY_AF_REF.loadSettingsFrom(settings);
        m_GTY_DP.loadSettingsFrom(settings);
        m_GTY_QUAL.loadSettingsFrom(settings);
        m_GTY_SEC_PL.loadSettingsFrom(settings);
        m_IGNORE_HOMO.loadSettingsFrom(settings);
        m_INFILE.loadSettingsFrom(settings);
        m_KGGSEQ.loadSettingsFrom(settings);
        m_OMIM_ANNO.loadSettingsFrom(settings);
        m_OUTFORMAT.loadSettingsFrom(settings);
        m_OUTPREFIX.loadSettingsFrom(settings);
        m_PEDFILE.loadSettingsFrom(settings);
        m_PUBMED.loadSettingsFrom(settings);
        m_SEQ_MQ.loadSettingsFrom(settings);
        m_SEQ_QUAL.loadSettingsFrom(settings);
        m_SEQ_SB.loadSettingsFrom(settings);
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


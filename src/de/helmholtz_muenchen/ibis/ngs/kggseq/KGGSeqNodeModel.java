package de.helmholtz_muenchen.ibis.ngs.kggseq;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
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
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.ngs.samsplitter.helpers.FileHelpers;

/**
 * This is the model implementation of KGGSeq.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class KGGSeqNodeModel extends HTExecutorNodeModel {
    
		
	/**
	 * Config Keys
	 */
	
//	Steps according to http://statgenpro.psychiatry.hku.hk/limx/kggseq/doc/Denovo.htm
	
//	Step 1
	public static final String CFGKEY_KGGSEQ_PATH 			= "KGGSEQ_PATH";
	public static final String CFGKEY_BUILDVER			  	= "BUILDVER";
	public static final String CFGKEY_INFILE			  	= "INFILE";
	public static final String CFGKEY_PEDFILE 			   	= "PEDFILE";
	public static final String CFGKEY_RESOURCE 			   	= "RESOURCE";
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
    private final SettingsModelString m_PEDFILE = new SettingsModelString(KGGSeqNodeModel.CFGKEY_PEDFILE, "");
    private final SettingsModelString m_RESOURCE = new SettingsModelString(KGGSeqNodeModel.CFGKEY_RESOURCE, "");
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
    private final SettingsModelOptionalString m_GENE_FEATURES = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_GENE_FEATURES, "0,1,2,3,4,5,6",true);
    private final SettingsModelOptionalString m_FILTER_COMMON = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_FILTER_COMMON,"hg19_1kg201204,hg19_dbsnp138,hg19_ESP6500AA,hg19_ESP6500EA", true);   
    private final SettingsModelBoolean m_DISEASE_CAUSING_PRED = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_DISEASE_CAUSING_PRED, false);
    private final SettingsModelBoolean m_OMIM_ANNO = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_OMIM_ANNO, false);
    private final SettingsModelBoolean m_CANDIDATE_PPI = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_CANDIDATE_PPI, false);
    private final SettingsModelString m_CANDIDATE_GENES = new SettingsModelString(KGGSeqNodeModel.CFGKEY_CANDIDATE_GENES, "");
    private final SettingsModelBoolean m_CANDIDATE_PATHWAYS = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_CANDIDATE_PATHWAYS, false);
    private final SettingsModelOptionalString m_PUBMED = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_PUBMED,"",false);

	
	//The Output Col Names
	public static final String OUT_COL1 = "Path2OutFile";
    
    /**
     * Constructor for the node model.
     */
    protected KGGSeqNodeModel() {
        super(1, 1);
        
        addSetting(m_BUILDVER);
        addSetting(m_CANDIDATE_GENES);
        addSetting(m_CANDIDATE_PATHWAYS);
        addSetting(m_CANDIDATE_PPI);
        addSetting(m_COMPOSITESUBJECTID);
        addSetting(m_DISEASE_CAUSING_PRED);
        addSetting(m_FILTER_COMMON);
        addSetting(m_GENE_FEATURES);
        addSetting(m_GENOTYPE_FILTER);
        addSetting(m_GTY_AF_ALT);
        addSetting(m_GTY_AF_HET);
        addSetting(m_GTY_AF_REF);
        addSetting(m_GTY_DP);
        addSetting(m_GTY_QUAL);
        addSetting(m_GTY_SEC_PL);
        addSetting(m_IGNORE_HOMO);
        addSetting(m_KGGSEQ);
        addSetting(m_OMIM_ANNO);
        addSetting(m_OUTFORMAT);
        addSetting(m_OUTPREFIX);
        addSetting(m_PEDFILE);
        addSetting(m_PUBMED);
        addSetting(m_SEQ_MQ);
        addSetting(m_SEQ_QUAL);
        addSetting(m_SEQ_SB);
        addSetting(m_RESOURCE);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String infile 		= inData[0].iterator().next().getCell(0).toString();
    	String BasePath 	= infile.substring(0,infile.lastIndexOf("/")+1);
    	String lockFile 	= BasePath+m_OUTPREFIX.getStringValue() + SuccessfulRunChecker.LOCK_ENDING;
    	String StdErrFile 	= BasePath+m_OUTPREFIX.getStringValue()+".stdErr.log";
    	String StdOutFile 	= BasePath+m_OUTPREFIX.getStringValue()+".stdOut.log";
    	String outfile 		= BasePath+m_OUTPREFIX.getStringValue()+".flt.xlsx";
    	
    	String[] command 		= createCommand(BasePath,infile);
    	boolean terminationState= false;
    	
    	//Check md5 sums to see if file has changed
    	String md5_oldOutFile 	= "";
    	String md5_newOutFile	= "";
    	
    	if(new File(outfile).exists()){
    		md5_oldOutFile = FileHelpers.getmd5Sum(outfile);
    		String lockCommand = "";
    		for (String s : command) {
    			lockCommand += s;
    		}
    		terminationState = SuccessfulRunChecker.hasTerminatedSuccessfully(new File(lockFile), lockCommand);
    	}
    	
    	
    	super.executeCommand(createCommand(BasePath,infile), exec, new File(lockFile),StdOutFile, StdErrFile);
    	
    	//Check if outfile was created--> execution was successful. 
    	if(!new File(outfile).exists()){
    		setWarningMessage("No outfile was created. Execution probably failed. Please check error logs. Removing .klock file.");
    		new File(lockFile).delete();
    	}else{
    		md5_newOutFile = FileHelpers.getmd5Sum(outfile);
    		
    		if(md5_oldOutFile.equals(md5_newOutFile) && !terminationState ){
        		setWarningMessage("Outfile already existed and did not change during execution. Execution probably failed! Please check error logs. Removing .klock file.");
        		new File(lockFile).delete();
    		}
    	}
    	
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(BasePath+m_OUTPREFIX.getStringValue())};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
		return new BufferedDataTable[]{outTable};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
    }

    
    /**
     * Creates the execution command of the KGGSeq Node
     * @return
     */
    private String[] createCommand(String BasePath, String infile){
    	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	command.add("java -jar");
    	command.add(m_KGGSEQ.getStringValue());
    	command.add("--buildver "+m_BUILDVER.getStringValue());
    	command.add("--resource "+m_RESOURCE.getStringValue());
    	command.add("--vcf-file "+infile);
    	command.add("--ped-file "+m_PEDFILE.getStringValue());
    	
    	if(m_COMPOSITESUBJECTID.getBooleanValue()){
    		command.add("--composite-subject-id");
    	}
    		
    	command.add("--out "+BasePath+m_OUTPREFIX.getStringValue());
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
    	
    	if(m_FILTER_COMMON.isActive()){
    		command.add("--db-filter "+m_FILTER_COMMON.getStringValue()+" --rare-allele-freq 0.05");
    	}
    	
    	if(m_DISEASE_CAUSING_PRED.getBooleanValue()){
    		command.add("--db-score dbnsfp --mendel-causing-predict all --filter-nondisease-variant");
    	}
    	
    	if(m_OMIM_ANNO.getBooleanValue()){
    		command.add("--genome-annot --omim-annot");
    	}
    	
    	if(m_CANDIDATE_GENES.isEnabled()){
    		if(!m_CANDIDATE_GENES.getStringValue().equals("")){
    			command.add("--candi-file "+m_CANDIDATE_GENES.getStringValue());
    		}
    		
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
    	command.add("--no-lib-check");
    	
    	
    	return new String[]{StringUtils.join(command, " ")};
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
    	
    	
		if(CompatibilityChecker.inputFileNotOk(m_KGGSEQ.getStringValue(), false)) {
			throw new InvalidSettingsException("Set path to KGGSeq jar!");
		}
    	
    	
    	if(CompatibilityChecker.getFirstIndexCellType(inSpecs[0], "VCFCell")!=0){
    		throw new InvalidSettingsException("Invalid input. No VCFCell in first column of input table.");
    	}
    	
        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec()})};
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


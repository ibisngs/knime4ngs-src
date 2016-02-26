package de.helmholtz_muenchen.ibis.ngs.vqsr;

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
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;

/**
 * This is the model implementation of VQSR.
 * 
 *
 * @author 
 */
public class VQSRNodeModel extends HTExecutorNodeModel {
    
	
	/**
	 * Logger
	 */
//	private static final NodeLogger LOGGER = NodeLogger.getLogger(VQSRNodeModel.class);
	
	/**
	 * CFGKEYS
	 */
	
    public static final String CFGKEY_GATK						= "gatk";
    public static final String CFGKEY_REF_GENOME				= "refgenome";
    public static final String CFGKEY_MODE						= "MODE";
    public static final String CFGKEY_AN						= "AN";
    public static final String CFGKEY_TRANCHE					= "Tranche";
    public static final String CFGKEY_GAUSS						= "Gaussians";
    public static final String CFGKEY_NT						= "NT";

    public static final String CFGKEY_RESOURCES_STRING_1		= "resources_1";
    public static final String CFGKEY_RESOURCES_STRING_2		= "resources_2";
    public static final String CFGKEY_RESOURCES_STRING_3		= "resources_3";
    public static final String CFGKEY_RESOURCES_STRING_4		= "resources_4";
    public static final String CFGKEY_RESOURCES_STRING_5		= "resources_5";

    public static final String CFGKEY_RESOURCES_BOOLEAN_1		= "resources_BOOLEAN_1";
    public static final String CFGKEY_RESOURCES_BOOLEAN_2		= "resources_BOOLEAN_2";
    public static final String CFGKEY_RESOURCES_BOOLEAN_3		= "resources_BOOLEAN_3";
    public static final String CFGKEY_RESOURCES_BOOLEAN_4		= "resources_BOOLEAN_4";
    public static final String CFGKEY_RESOURCES_BOOLEAN_5		= "resources_BOOLEAN_5";

    public static final String CFGKEY_RESOURCES_FILE_1			= "resources_FILE_1";
    public static final String CFGKEY_RESOURCES_FILE_2			= "resources_FILE_2";
    public static final String CFGKEY_RESOURCES_FILE_3			= "resources_FILE_3";
    public static final String CFGKEY_RESOURCES_FILE_4			= "resources_FILE_4";
    public static final String CFGKEY_RESOURCES_FILE_5			= "resources_FILE_5";

    public static final String CFGKEY_TS_FILTER					= "TS_FILTER";
    
    public static final String CFGKEY_OPT_VAR_RECAL 			= "var_recal_opts";
    public static final String CFGKEY_OPT_APPLY_RECAL			= "apply_recal_opts";
    
    public static final String DEFAULT_MODE						= "SNP";
    public static final String DEFAULT_SNP_AN					= "-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff -an MQ";
    public static final String DEFAULT_INDEL_AN					= "-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff";
    public static final String DEFAULT_TRANCHES					= "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0";
    public static final double DEFAULT_TS_FILTER_SNP			= 99.5;
    public static final double DEFAULT_TS_FILTER_INDEL			= 99.0;
    
	/**
	 * Models
	 */
    
    private final SettingsModelString m_GATK = new SettingsModelString(VQSRNodeModel.CFGKEY_GATK, "");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(VQSRNodeModel.CFGKEY_REF_GENOME, "");
    private final SettingsModelString m_MODE = new SettingsModelString(VQSRNodeModel.CFGKEY_MODE, VQSRNodeModel.DEFAULT_MODE);
    private final SettingsModelString m_TRANCHE = new SettingsModelString(VQSRNodeModel.CFGKEY_TRANCHE, VQSRNodeModel.DEFAULT_TRANCHES);
    private final SettingsModelString m_AN = new SettingsModelString(VQSRNodeModel.CFGKEY_AN, VQSRNodeModel.DEFAULT_SNP_AN);
    private final SettingsModelIntegerBounded m_GAUSSIANS = new SettingsModelIntegerBounded(VQSRNodeModel.CFGKEY_GAUSS,8,1,Integer.MAX_VALUE);
    private final SettingsModelIntegerBounded m_NT = new SettingsModelIntegerBounded(VQSRNodeModel.CFGKEY_NT,1,1,Integer.MAX_VALUE);

    private final SettingsModelString m_RESOURCES_STRING_1 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_1, "hapmap,known=false,training=true,truth=true,prior=15.0");
    private final SettingsModelString m_RESOURCES_STRING_2 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_2, "omni,known=false,training=true,truth=true,prior=12.0");
    private final SettingsModelString m_RESOURCES_STRING_3 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_3, "1000G,known=false,training=true,truth=false,prior=10.0");
    private final SettingsModelString m_RESOURCES_STRING_4 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_4, "dbsnp,known=true,training=false,truth=false,prior=2.0");
    private final SettingsModelString m_RESOURCES_STRING_5 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_5, "mills,known=false,training=true,truth=true,prior=12.0");
    
    //TODO remove for publication
    private final SettingsModelString m_RESOURCES_FILE_1 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_1, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/hapmap_3.3.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_2 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_2, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/1000G_omni2.5.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_3 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_3, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/1000G_phase1.snps.high_confidence.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_4 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_4, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/dbsnp_138.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_5 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_5, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/Mills_and_1000G_gold_standard.indels.hg19.vcf");

    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_1 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_1, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_2 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_2, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_3 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_3, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_4 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_4, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_5 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_5, false);

    private final SettingsModelDoubleBounded m_TS_FILTER = new SettingsModelDoubleBounded(VQSRNodeModel.CFGKEY_TS_FILTER,VQSRNodeModel.DEFAULT_TS_FILTER_SNP,1,100);

    private final SettingsModelOptionalString m_OPT_VAR_RECAL = new SettingsModelOptionalString(VQSRNodeModel.CFGKEY_OPT_VAR_RECAL,"",false);
    private final SettingsModelOptionalString m_OPT_APPLY_RECAL = new SettingsModelOptionalString(VQSRNodeModel.CFGKEY_OPT_APPLY_RECAL,"",false);
    
	//The Output Col Names
	public static final String OUT_COL1 = "VQSR VARIANTS";
    
    /**
     * Constructor for the node model.
     */
    protected VQSRNodeModel() {
    
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
    	//Check if path to GATK is available
    	String PATH2GATK = "";
    	int GATK_INDEX = inData[0].getDataTableSpec().findColumnIndex("Path2GATKFile");
    	if(GATK_INDEX!=-1){
    		PATH2GATK = inData[0].iterator().next().getCell(GATK_INDEX).toString();
    	}
    	//Check if path to GATK is available
    	String PATH2REFSEQ = "";
    	int REFSEQ_INDEX = inData[0].getDataTableSpec().findColumnIndex("Path2SEQFile");
    	if(REFSEQ_INDEX!=-1){
    		PATH2REFSEQ = inData[0].iterator().next().getCell(REFSEQ_INDEX).toString();
    	}
    	
    	//Check Infile
    	String INFILE;
    	try{
    		INFILE = inData[0].iterator().next().getCell(0).toString();
    		if(!INFILE.endsWith(".vcf")){
    			throw new InvalidSettingsException("First Cell of input table has to be the path to VCF Infile but it is "+INFILE);
    		}
    	}catch(IndexOutOfBoundsException e){
    			throw new InvalidSettingsException("First Cell of input table has to be the path to VCF Infile but it is empty.");
    	}
    	
    	String outFile;
    	if(INFILE.endsWith(".vcf")){
    		outFile = INFILE.replace(".vcf","_"+m_MODE.getStringValue()+"_VQSR.vcf");
    	}else{
    		outFile = INFILE+"_"+m_MODE.getStringValue()+"_VQSR.vcf";
    	}
    	
    	File lockFile1 = new File(IO.replaceFileExtension(outFile,"var_recal" + SuccessfulRunChecker.LOCK_ENDING));
    	File lockFile2 = new File(IO.replaceFileExtension(outFile,"apply_recal" + SuccessfulRunChecker.LOCK_ENDING));
    	String stdOut1 = IO.replaceFileExtension(outFile,"var_recal.stdOut");
    	String stdErr1 = IO.replaceFileExtension(outFile,"var_recal.stdErr");
    	String stdOut2 = IO.replaceFileExtension(outFile,"apply_recal.stdOut");
    	String stdErr2 = IO.replaceFileExtension(outFile,"apply_recal.stdErr");
    	
    	//Execute
//    	Executor.executeCommand(createRecalibrationCommand(PATH2GATK,PATH2REFSEQ,INFILE),exec,LOGGER,stdOut,stdErr);
    	super.executeCommand(createRecalibrationCommand(PATH2GATK,PATH2REFSEQ,INFILE),exec,lockFile1,stdOut1,stdErr1);

//    	Executor.executeCommand(applyRecalibrationCommand(PATH2GATK,PATH2REFSEQ,INFILE),exec,LOGGER,stdOut,stdErr);
    	super.executeCommand(applyRecalibrationCommand(PATH2GATK,PATH2REFSEQ,INFILE),exec,lockFile2,stdOut2,stdErr2);

    	
    	/**
    	 * OUTPUT
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(outFile)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	
        return new BufferedDataTable[]{outTable};
    }

    /**
     *Creates Recalibration Command
     * @return
     */
    private String[] createRecalibrationCommand(String PATH2GATK, String PATH2REFSEQ, String INFILE){
    	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	String GATK;
    	if(!PATH2GATK.equals("")){
    		GATK = PATH2GATK;
    	}else{
    		GATK = m_GATK.getStringValue();
    	}
    	
    	String REFSEQ;
    	if(!PATH2REFSEQ.equals("")){
    		REFSEQ = PATH2REFSEQ;
    	}else{
    		REFSEQ = m_REF_GENOME.getStringValue();
    	}
    
    	String infile		= INFILE;
    	String recalFile 	= "";
    	String tranchesFile	= "";
    	String plotFile		= "";
    	if(infile.endsWith(".vcf")){
    		recalFile = infile.replace(".vcf","_"+m_MODE.getStringValue()+ "_VQSR.recal");
    		tranchesFile = infile.replace(".vcf","_"+m_MODE.getStringValue() +"_VQSR.tranches");
    		plotFile = infile.replace(".vcf","_"+m_MODE.getStringValue()+"_VQSR.plots.R");
    	}else{
    		recalFile = infile+"_"+m_MODE.getStringValue()+"_VQSR.recal";
    		tranchesFile = infile+"_"+m_MODE.getStringValue()+"_VQSR.tranches";
    		plotFile = infile+"_"+m_MODE.getStringValue()+"_VQSR.plots.R";
    	}
    	
    	command.add("java -jar");
    	command.add(GATK);
    	command.add("-T VariantRecalibrator");
    	command.add("-R "+REFSEQ);
    	command.add("-input "+infile);
    	
    	//Add resources
    	if(m_RESOURCES_BOOLEAN_1.getBooleanValue()){
    		command.add("-resource:"+m_RESOURCES_STRING_1.getStringValue()+" "+m_RESOURCES_FILE_1.getStringValue());
    	}
    	if(m_RESOURCES_BOOLEAN_2.getBooleanValue()){
    		command.add("-resource:"+m_RESOURCES_STRING_2.getStringValue()+" "+m_RESOURCES_FILE_2.getStringValue());
    	}
    	if(m_RESOURCES_BOOLEAN_3.getBooleanValue()){
    		command.add("-resource:"+m_RESOURCES_STRING_3.getStringValue()+" "+m_RESOURCES_FILE_3.getStringValue());
    	}
    	if(m_RESOURCES_BOOLEAN_4.getBooleanValue()){
    		command.add("-resource:"+m_RESOURCES_STRING_4.getStringValue()+" "+m_RESOURCES_FILE_4.getStringValue());
    	}
    	
    	if(m_RESOURCES_BOOLEAN_5.getBooleanValue()){
    		command.add("-resource:"+m_RESOURCES_STRING_5.getStringValue()+" "+m_RESOURCES_FILE_5.getStringValue());
    	}
    	
    	command.add(m_AN.getStringValue());
    	
    	command.add("-mode "+m_MODE.getStringValue());
    	
    	command.add("-recalFile "+recalFile);
    	command.add("-tranchesFile "+tranchesFile);
    	command.add("-rscriptFile "+plotFile);
    	command.add("-nt "+m_NT.getIntValue());
    	command.add(m_TRANCHE.getStringValue());
    	command.add("--maxGaussians "+m_GAUSSIANS.getIntValue());
    	command.add(m_OPT_VAR_RECAL.getStringValue());
    	
    	return new String[]{StringUtils.join(command, " ")};
    }
    
    /**
     *Creates Recalibration Command
     * @return
     */
    private String[] applyRecalibrationCommand(String PATH2GATK, String PATH2REFSEQ,String INFILE){
    	
    	ArrayList<String> command = new ArrayList<String>();
    	
    	String GATK;
    	if(!PATH2GATK.equals("")){
    		GATK = PATH2GATK;
    	}else{
    		GATK = m_GATK.getStringValue();
    	}
    	String REFSEQ;
    	if(!PATH2REFSEQ.equals("")){
    		REFSEQ = PATH2REFSEQ;
    	}else{
    		REFSEQ = m_REF_GENOME.getStringValue();
    	}
    	
    	//Process Infile Name to obtain output file names
    	String infile 		= INFILE;
    	String recalFile 	= "";
    	String tranchesFile	= "";
    	String outFile		= "";
    	if(infile.endsWith(".vcf")){
    		recalFile = infile.replace(".vcf","_"+m_MODE.getStringValue()+"_VQSR.recal");
    		tranchesFile = infile.replace(".vcf","_"+m_MODE.getStringValue()+"_VQSR.tranches");
    		outFile = infile.replace(".vcf","_"+m_MODE.getStringValue()+"_VQSR.vcf");
    	}else{
    		recalFile = infile+"_"+m_MODE.getStringValue()+"_VQSR.recal";
    		tranchesFile = infile+"_"+m_MODE.getStringValue()+"_VQSR.tranches";
    		outFile = infile+"_"+m_MODE.getStringValue()+"_VQSR.vcf";
    	}
    	
    	command.add("java -jar");
    	command.add(GATK);
    	command.add("-T ApplyRecalibration");
    	command.add("-R "+REFSEQ);
    	command.add("-input "+infile); 	
       	command.add("-mode "+m_MODE.getStringValue());
    	
    	command.add("-recalFile "+recalFile);
    	command.add("-tranchesFile "+tranchesFile);
    	command.add("-o "+outFile);

    	command.add("-ts_filter_level "+m_TS_FILTER.getDoubleValue());
    	command.add(m_OPT_APPLY_RECAL.getStringValue());
    	    	
    	return new String[]{StringUtils.join(command, " ")};
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

    	
//    	if(inSpecs[0].containsName("Path2GATKFile")){
//    		m_GATK.setEnabled(false);
//    		m_GATK.setStringValue("GATK Path already set");
//    		LOGGER.info("Found GATK Path in inData");
//    	}else{
//    		LOGGER.warn("No path to GATK found in inData. User input required.");
//    	}
//    	
//    	if(inSpecs[0].containsName("Path2SEQFile")){
//    		m_REF_GENOME.setEnabled(false);
//    		m_REF_GENOME.setStringValue("Reference Sequence already set");
//    		LOGGER.info("Found Reference Sequence in inData");
//    	}else{
//    		LOGGER.warn("No path to Reference Sequence found in inData. User input required.");
//    	}
    	

    	return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, VCFCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	super.saveSettingsTo(settings);
         m_AN.saveSettingsTo(settings);
         m_GATK.saveSettingsTo(settings);
         m_MODE.saveSettingsTo(settings);
         m_NT.saveSettingsTo(settings);
         m_GAUSSIANS.saveSettingsTo(settings);
         m_REF_GENOME.saveSettingsTo(settings);
         m_TRANCHE.saveSettingsTo(settings);
         
         m_RESOURCES_BOOLEAN_1.saveSettingsTo(settings);
         m_RESOURCES_BOOLEAN_2.saveSettingsTo(settings);
         m_RESOURCES_BOOLEAN_3.saveSettingsTo(settings);
         m_RESOURCES_BOOLEAN_4.saveSettingsTo(settings);
         m_RESOURCES_BOOLEAN_5.saveSettingsTo(settings);
         m_RESOURCES_FILE_1.saveSettingsTo(settings);
         m_RESOURCES_FILE_2.saveSettingsTo(settings);
         m_RESOURCES_FILE_3.saveSettingsTo(settings);
         m_RESOURCES_FILE_4.saveSettingsTo(settings);
         m_RESOURCES_FILE_5.saveSettingsTo(settings);
         m_RESOURCES_STRING_1.saveSettingsTo(settings);
         m_RESOURCES_STRING_2.saveSettingsTo(settings);
         m_RESOURCES_STRING_3.saveSettingsTo(settings);
         m_RESOURCES_STRING_4.saveSettingsTo(settings);
         m_RESOURCES_STRING_5.saveSettingsTo(settings);
         m_TS_FILTER.saveSettingsTo(settings);
         
         m_OPT_VAR_RECAL.saveSettingsTo(settings);
         m_OPT_APPLY_RECAL.saveSettingsTo(settings);
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
        m_AN.loadSettingsFrom(settings);
        m_GATK.loadSettingsFrom(settings);
        m_MODE.loadSettingsFrom(settings);
        m_NT.loadSettingsFrom(settings);
        m_GAUSSIANS.loadSettingsFrom(settings);
        m_REF_GENOME.loadSettingsFrom(settings);
        m_TRANCHE.loadSettingsFrom(settings);
        
        m_RESOURCES_BOOLEAN_1.loadSettingsFrom(settings);
        m_RESOURCES_BOOLEAN_2.loadSettingsFrom(settings);
        m_RESOURCES_BOOLEAN_3.loadSettingsFrom(settings);
        m_RESOURCES_BOOLEAN_4.loadSettingsFrom(settings);
        m_RESOURCES_BOOLEAN_5.loadSettingsFrom(settings);
        m_RESOURCES_FILE_1.loadSettingsFrom(settings);
        m_RESOURCES_FILE_2.loadSettingsFrom(settings);
        m_RESOURCES_FILE_3.loadSettingsFrom(settings);
        m_RESOURCES_FILE_4.loadSettingsFrom(settings);
        m_RESOURCES_FILE_5.loadSettingsFrom(settings);
        m_RESOURCES_STRING_1.loadSettingsFrom(settings);
        m_RESOURCES_STRING_2.loadSettingsFrom(settings);
        m_RESOURCES_STRING_3.loadSettingsFrom(settings);
        m_RESOURCES_STRING_4.loadSettingsFrom(settings);
        m_RESOURCES_STRING_5.loadSettingsFrom(settings);
        m_TS_FILTER.loadSettingsFrom(settings);
        
        m_OPT_VAR_RECAL.loadSettingsFrom(settings);
        m_OPT_APPLY_RECAL.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	super.validateSettings(settings);
        m_AN.validateSettings(settings);
        m_GATK.validateSettings(settings);
        m_MODE.validateSettings(settings);
        m_NT.validateSettings(settings);
        m_GAUSSIANS.validateSettings(settings);
        m_REF_GENOME.validateSettings(settings);
        m_TRANCHE.validateSettings(settings);
        
        m_RESOURCES_BOOLEAN_1.validateSettings(settings);
        m_RESOURCES_BOOLEAN_2.validateSettings(settings);
        m_RESOURCES_BOOLEAN_3.validateSettings(settings);
        m_RESOURCES_BOOLEAN_4.validateSettings(settings);
        m_RESOURCES_BOOLEAN_5.validateSettings(settings);
        m_RESOURCES_FILE_1.validateSettings(settings);
        m_RESOURCES_FILE_2.validateSettings(settings);
        m_RESOURCES_FILE_3.validateSettings(settings);
        m_RESOURCES_FILE_4.validateSettings(settings);
        m_RESOURCES_FILE_5.validateSettings(settings);
        m_RESOURCES_STRING_1.validateSettings(settings);
        m_RESOURCES_STRING_2.validateSettings(settings);
        m_RESOURCES_STRING_3.validateSettings(settings);
        m_RESOURCES_STRING_4.validateSettings(settings);
        m_RESOURCES_STRING_5.validateSettings(settings);
        m_TS_FILTER.validateSettings(settings);
        
        m_OPT_VAR_RECAL.validateSettings(settings);
        m_OPT_APPLY_RECAL.validateSettings(settings);
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


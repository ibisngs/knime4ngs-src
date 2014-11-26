package de.helmholtz_muenchen.ibis.ngs.vqsr;

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
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of VQSR.
 * 
 *
 * @author 
 */
public class VQSRNodeModel extends NodeModel {
    
	
	/**
	 * Logger
	 */
	private static final NodeLogger LOGGER = NodeLogger.getLogger(VQSRNodeModel.class);
	
	/**
	 * CFGKEYS
	 */
	
    public static final String CFGKEY_GATK						= "gatk";
    public static final String CFGKEY_REF_GENOME				= "refgenome";
    public static final String CFGKEY_INFILE					= "infile";
    public static final String CFGKEY_MODE						= "MODE";
    public static final String CFGKEY_AN						= "AN";
    public static final String CFGKEY_TRANCHE					= "Tranche";
    public static final String CFGKEY_NT						= "NT";

    public static final String CFGKEY_RESOURCES_STRING_1		= "resources_1";
    public static final String CFGKEY_RESOURCES_STRING_2		= "resources_2";
    public static final String CFGKEY_RESOURCES_STRING_3		= "resources_3";
    public static final String CFGKEY_RESOURCES_STRING_4		= "resources_4";
    public static final String CFGKEY_RESOURCES_BOOLEAN_1		= "resources_BOOLEAN_1";
    public static final String CFGKEY_RESOURCES_BOOLEAN_2		= "resources_BOOLEAN_2";
    public static final String CFGKEY_RESOURCES_BOOLEAN_3		= "resources_BOOLEAN_3";
    public static final String CFGKEY_RESOURCES_BOOLEAN_4		= "resources_BOOLEAN_4";
    public static final String CFGKEY_RESOURCES_FILE_1			= "resources_FILE_1";
    public static final String CFGKEY_RESOURCES_FILE_2			= "resources_FILE_2";
    public static final String CFGKEY_RESOURCES_FILE_3			= "resources_FILE_3";
    public static final String CFGKEY_RESOURCES_FILE_4			= "resources_FILE_4";
    
    public static final String CFGKEY_TS_FILTER					= "TS_FILTER";
    
	/**
	 * Models
	 */
    
    private final SettingsModelString m_GATK = new SettingsModelString(VQSRNodeModel.CFGKEY_GATK, "");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(VQSRNodeModel.CFGKEY_REF_GENOME, "");
//    private final SettingsModelString m_INFILE = new SettingsModelString(VQSRNodeModel.CFGKEY_INFILE, "");
    private final SettingsModelString m_MODE = new SettingsModelString(VQSRNodeModel.CFGKEY_MODE, "SNP");
    private final SettingsModelString m_TRANCHE = new SettingsModelString(VQSRNodeModel.CFGKEY_TRANCHE, "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0");
    private final SettingsModelString m_AN = new SettingsModelString(VQSRNodeModel.CFGKEY_AN, "-an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum");
    private final SettingsModelIntegerBounded m_NT = new SettingsModelIntegerBounded(VQSRNodeModel.CFGKEY_NT,1,1,Integer.MAX_VALUE);

    private final SettingsModelString m_RESOURCES_STRING_1 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_1, "hapmap,known=false,training=true,truth=true,prior=15.0");
    private final SettingsModelString m_RESOURCES_STRING_2 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_2, "omni,known=false,training=true,truth=true,prior=12.0");
    private final SettingsModelString m_RESOURCES_STRING_3 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_3, "1000G,known=false,training=true,truth=false,prior=10.0");
    private final SettingsModelString m_RESOURCES_STRING_4 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_4, "dbsnp,known=true,training=false,truth=false,prior=2.0");
    
    private final SettingsModelString m_RESOURCES_FILE_1 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_1, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/hapmap_3.3.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_2 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_2, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/1000G_omni2.5.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_3 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_3, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/1000G_phase1.snps.high_confidence.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_4 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_4, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/dbsnp_138.hg19.vcf");
    
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_1 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_1, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_2 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_2, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_3 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_3, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_4 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_4, true);

    private final SettingsModelDoubleBounded m_TS_FILTER = new SettingsModelDoubleBounded(VQSRNodeModel.CFGKEY_TS_FILTER,99.9,1,100);

    /**
     * Constructor for the node model.
     */
    protected VQSRNodeModel() {
    
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
    	
    	String INFILE = inData[0].iterator().next().getCell(0).toString();
    	
    	Executor.executeCommand(createRecalibrationCommand(PATH2GATK,PATH2REFSEQ,INFILE),exec,LOGGER,stdOut,stdErr);
    	Executor.executeCommand(applyRecalibrationCommand(PATH2GATK,PATH2REFSEQ,INFILE),exec,LOGGER,stdOut,stdErr);
    	
        // TODO: Return a BufferedDataTable for each output port 
        return inData;
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
    	
    	//Process Infile Name to obtain output file names
//    	String infile 		= m_INFILE.getStringValue();
    	String infile		= INFILE;
    	
    	
    	String recalFile 	= "";
    	String tranchesFile	= "";
    	String plotFile		= "";
    	if(infile.endsWith(".vcf")){
    		recalFile = infile.replace(".vcf","_VQSR.recal");
    		tranchesFile = infile.replace(".vcf","_VQSR.tranches");
    		plotFile = infile.replace(".vcf","_VQSR.plots.R");
    	}else{
    		recalFile = infile+"_VQSR.recal";
    		tranchesFile = infile+"_VQSR.tranches";
    		plotFile = infile+"_VQSR.plots.R";
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
    	
    	command.add(m_AN.getStringValue());
    	
    	command.add("-mode "+m_MODE.getStringValue());
    	
    	command.add("-recalFile "+recalFile);
    	command.add("-tranchesFile "+tranchesFile);
    	command.add("-rscriptFile "+plotFile);
    	command.add(m_TRANCHE.getStringValue());
    	
    	System.out.println(StringUtils.join(command, " "));
    	
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
    		recalFile = infile.replace(".vcf","_VQSR.recal");
    		tranchesFile = infile.replace(".vcf","_VQSR.tranches");
    		outFile = infile.replace(".vcf","_"+m_MODE.getStringValue()+"_VQSR.vcf");
    	}else{
    		recalFile = infile+"_VQSR.recal";
    		tranchesFile = infile+"_VQSR.tranches";
    		outFile = infile+"_VQSR.vcf";
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
    	
    	System.out.println(StringUtils.join(command, " "));
    	
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

    	if(inSpecs[0].containsName("Path2GATKFile")){
    		m_GATK.setEnabled(false);
    		m_GATK.setStringValue("GATK Path already set");
    		System.out.println("Found GATK Path in inData");
    	}else{
    		System.out.println("No path to GATK found in inData. User input required.");
    	}
    	
    	if(inSpecs[0].containsName("Path2SEQFile")){
    		m_REF_GENOME.setEnabled(false);
    		m_REF_GENOME.setStringValue("Reference Sequence already set");
    		System.out.println("Found Reference Sequence in inData");
    	}else{
    		System.out.println("No path to Reference Sequence found in inData. User input required.");
    	}
    	
    	
        // TODO: generated method stub
        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         m_AN.saveSettingsTo(settings);
         m_GATK.saveSettingsTo(settings);
//         m_INFILE.saveSettingsTo(settings);
         m_MODE.saveSettingsTo(settings);
         m_NT.saveSettingsTo(settings);
         m_REF_GENOME.saveSettingsTo(settings);
         m_TRANCHE.saveSettingsTo(settings);
         
         m_RESOURCES_BOOLEAN_1.saveSettingsTo(settings);
         m_RESOURCES_BOOLEAN_2.saveSettingsTo(settings);
         m_RESOURCES_BOOLEAN_3.saveSettingsTo(settings);
         m_RESOURCES_BOOLEAN_4.saveSettingsTo(settings);
         m_RESOURCES_FILE_1.saveSettingsTo(settings);
         m_RESOURCES_FILE_2.saveSettingsTo(settings);
         m_RESOURCES_FILE_3.saveSettingsTo(settings);
         m_RESOURCES_FILE_4.saveSettingsTo(settings);
         m_RESOURCES_STRING_1.saveSettingsTo(settings);
         m_RESOURCES_STRING_2.saveSettingsTo(settings);
         m_RESOURCES_STRING_3.saveSettingsTo(settings);
         m_RESOURCES_STRING_4.saveSettingsTo(settings);
         m_TS_FILTER.saveSettingsTo(settings);
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_AN.loadSettingsFrom(settings);
        m_GATK.loadSettingsFrom(settings);
//        m_INFILE.loadSettingsFrom(settings);
        m_MODE.loadSettingsFrom(settings);
        m_NT.loadSettingsFrom(settings);
        m_REF_GENOME.loadSettingsFrom(settings);
        m_TRANCHE.loadSettingsFrom(settings);
        
        m_RESOURCES_BOOLEAN_1.loadSettingsFrom(settings);
        m_RESOURCES_BOOLEAN_2.loadSettingsFrom(settings);
        m_RESOURCES_BOOLEAN_3.loadSettingsFrom(settings);
        m_RESOURCES_BOOLEAN_4.loadSettingsFrom(settings);
        m_RESOURCES_FILE_1.loadSettingsFrom(settings);
        m_RESOURCES_FILE_2.loadSettingsFrom(settings);
        m_RESOURCES_FILE_3.loadSettingsFrom(settings);
        m_RESOURCES_FILE_4.loadSettingsFrom(settings);
        m_RESOURCES_STRING_1.loadSettingsFrom(settings);
        m_RESOURCES_STRING_2.loadSettingsFrom(settings);
        m_RESOURCES_STRING_3.loadSettingsFrom(settings);
        m_RESOURCES_STRING_4.loadSettingsFrom(settings);
        m_TS_FILTER.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_AN.validateSettings(settings);
        m_GATK.validateSettings(settings);
//        m_INFILE.validateSettings(settings);
        m_MODE.validateSettings(settings);
        m_NT.validateSettings(settings);
        m_REF_GENOME.validateSettings(settings);
        m_TRANCHE.validateSettings(settings);
        
        m_RESOURCES_BOOLEAN_1.validateSettings(settings);
        m_RESOURCES_BOOLEAN_2.validateSettings(settings);
        m_RESOURCES_BOOLEAN_3.validateSettings(settings);
        m_RESOURCES_BOOLEAN_4.validateSettings(settings);
        m_RESOURCES_FILE_1.validateSettings(settings);
        m_RESOURCES_FILE_2.validateSettings(settings);
        m_RESOURCES_FILE_3.validateSettings(settings);
        m_RESOURCES_FILE_4.validateSettings(settings);
        m_RESOURCES_STRING_1.validateSettings(settings);
        m_RESOURCES_STRING_2.validateSettings(settings);
        m_RESOURCES_STRING_3.validateSettings(settings);
        m_RESOURCES_STRING_4.validateSettings(settings);
        m_TS_FILTER.validateSettings(settings);
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


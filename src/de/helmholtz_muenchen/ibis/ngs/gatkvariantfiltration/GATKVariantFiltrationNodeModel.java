package de.helmholtz_muenchen.ibis.ngs.gatkvariantfiltration;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

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


/**
 * This is the model implementation of GATKVariantFiltration.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKVariantFiltrationNodeModel extends NodeModel {
    
	
	//--filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0
	
	
	/**
	 * Config Keys
	 */
	public static final String CFGKEY_GATK_PATH = "GATK_PATH";
	public static final String CFGKEY_INFILE = "INFILE";
	public static final String CFGKEY_REF_GENOME = "REFGENOME";
	public static final String CFGKEY_QUAL = "Quality score";
	public static final String CFGKEY_DP = "depth of reads ";
	public static final String CFGKEY_AD = "allele depth ";
	public static final String CFGKEY_QD = "QualByDepth";
	public static final String CFGKEY_FS = "Fisherstrand";
	public static final String CFGKEY_MQ = "MappingQuality";
	public static final String CFGKEY_HS = "HaplotypeScore";
	public static final String CFGKEY_MQR = "MappingQualityRankSum";
	public static final String CFGKEY_RPR = "ReadPosRankSum";
	public static final String CFGKEY_FilterString = "FilterString";
	public static final String CFGKEY_FilterName = "FilterName";
	
	
	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_GATK = new SettingsModelString(CFGKEY_GATK_PATH, "");
    private final SettingsModelString m_INFILE = new SettingsModelString(CFGKEY_INFILE, "");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(CFGKEY_REF_GENOME, "");
    
	private final SettingsModelString m_QUAL= new SettingsModelString(CFGKEY_QUAL, "50.0");
	private final SettingsModelString m_DP= new SettingsModelString(CFGKEY_DP, "20.0");
	private final SettingsModelString m_AD= new SettingsModelString(CFGKEY_AD, "5.0");
	private final SettingsModelString m_QD= new SettingsModelString(CFGKEY_QD, "2.0");
	private final SettingsModelString m_FS= new SettingsModelString(CFGKEY_FS, "60.0");
	private final SettingsModelString m_MQ= new SettingsModelString(CFGKEY_MQ, "40.0");
	private final SettingsModelString m_HS= new SettingsModelString(CFGKEY_HS, "13.0");
	private final SettingsModelString m_MQR= new SettingsModelString(CFGKEY_MQR, "-12.5");
	private final SettingsModelString m_RPR= new SettingsModelString(CFGKEY_RPR, "-8.0");
	private final SettingsModelString m_FilterString = new SettingsModelString(CFGKEY_FilterString, "-");
	private final SettingsModelString m_FilterName = new SettingsModelString(CFGKEY_FilterName, "-");
	
	
	/**
	 * Logger
	 */
	private static final NodeLogger LOGGER = NodeLogger.getLogger(GATKVariantFiltrationNodeModel.class);
	
	/**
	 * String that holds the complete filter options
	 */
	private static StringBuffer filterString = new StringBuffer();
	
    /**
     * Constructor for the node model.
     */
    protected GATKVariantFiltrationNodeModel() {
    
    	
    	
    	
        // TODO: Specify the amount of input and output ports needed.
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	ArrayList<String> command = new ArrayList<String>();
 //$com = "java -Xmx4g -jar ".$TOOL_PATH.$GATK_VERSION."/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_GENOME -V $INFILE
 //   	-o $outfile  --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0
 //   	|| MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0\" --filterName \"GATKFilter\"";
    	
    	command.add("java");
    	command.add("-Xmx4g -jar "+m_GATK.getStringValue());
    	command.add("-T VariantFiltration");
    	
    	command.add("-R "+m_REF_GENOME.getStringValue());
    	command.add("-V "+m_INFILE.getStringValue());
    	
    	String OUTFILE = m_INFILE.getStringValue().replaceAll(".vcf", "GATKfiltered.vcf");
    	command.add("-o "+OUTFILE);
    	command.add("--filterExpression");
  
    	/**
    	 * Create Filter String
    	 */
    	addToFilterString(m_QUAL);
    	addToFilterString(m_DP);
    	addToFilterString(m_AD);
    	addToFilterString(m_QD);
    	addToFilterString(m_FS);
    	addToFilterString(m_MQ);
    	addToFilterString(m_HS);
    	addToFilterString(m_MQR);
    	addToFilterString(m_RPR);
    	addToFilterString(m_FilterString);

    	
    	command.add(filterString.toString());
    	
    	if(m_FilterName.isEnabled()){
    		command.add("--filterName");
    		command.add(m_FilterName.getStringValue());
    	}else{
    		command.add("--filterName");
    		command.add("GATKVariantFiltration");
    	}
    	

        return new BufferedDataTable[]{};
    }

    /**
     * Adds the values of SettingsModelString to the filterString if model is enabled
     * @param toAdd
     */
    private void addToFilterString(SettingsModelString toAdd){
    	if(toAdd.isEnabled()){
        	if(filterString.length()==0){
        		filterString.append(toAdd.getStringValue());
        	}else{
        		filterString.append(" || "+toAdd.getStringValue());
        	}
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

        // TODO: generated method stub
        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	 m_GATK.saveSettingsTo(settings);
    	 m_INFILE.saveSettingsTo(settings);
    	 m_REF_GENOME.saveSettingsTo(settings);
         m_AD.saveSettingsTo(settings);
         m_DP.saveSettingsTo(settings);
         m_FS.saveSettingsTo(settings);
         m_HS.saveSettingsTo(settings);
         m_MQ.saveSettingsTo(settings);
         m_MQR.saveSettingsTo(settings);
         m_QD.saveSettingsTo(settings);
         m_QUAL.saveSettingsTo(settings);
         m_RPR.saveSettingsTo(settings);
         m_FilterString.saveSettingsTo(settings);
         m_FilterName.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_GATK.loadSettingsFrom(settings);
   	 	m_INFILE.loadSettingsFrom(settings);
   	 	m_REF_GENOME.loadSettingsFrom(settings);
    	m_AD.loadSettingsFrom(settings);
        m_DP.loadSettingsFrom(settings);
        m_FS.loadSettingsFrom(settings);
        m_HS.loadSettingsFrom(settings);
        m_MQ.loadSettingsFrom(settings);
        m_MQR.loadSettingsFrom(settings);
        m_QD.loadSettingsFrom(settings);
        m_QUAL.loadSettingsFrom(settings);
        m_RPR.loadSettingsFrom(settings);
        m_FilterString.loadSettingsFrom(settings);
        m_FilterName.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	m_GATK.validateSettings(settings);
    	m_INFILE.validateSettings(settings);
    	m_REF_GENOME.validateSettings(settings);
    	m_AD.validateSettings(settings);
        m_DP.validateSettings(settings);
        m_FS.validateSettings(settings);
        m_HS.validateSettings(settings);
        m_MQ.validateSettings(settings);
        m_MQR.validateSettings(settings);
        m_QD.validateSettings(settings);
        m_QUAL.validateSettings(settings);
        m_RPR.validateSettings(settings);
        m_FilterString.validateSettings(settings);
        m_FilterName.validateSettings(settings);
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


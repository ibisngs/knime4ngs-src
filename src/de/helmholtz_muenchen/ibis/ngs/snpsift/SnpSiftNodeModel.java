package de.helmholtz_muenchen.ibis.ngs.snpsift;

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

import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;



/**
 * This is the model implementation of SnpSift.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class SnpSiftNodeModel extends NodeModel {
    
	public static final String CFGKEY_SNPEFF_FOLDER = "snpeff_folder";
	public static final String CFGKEY_INVCF = "invcf";
	public static final String CFGKEY_METHOD = "method";
	
	//Filter
	public static final String CFGKEY_FILTERSTRING = "filterstring";
	public static final String CFGKEY_FILTERQUAL = "filterqual";
	public static final String CFGKEY_FILTERCOVERAGE = "filtercoverage";
	public static final String CFGKEY_FILTERQUALBOOL = "filterqualbool";
	public static final String CFGKEY_FILTERCOVERAGEBOOL = "filtercoveragebool";
	
	//Annotate
	public static final String CFGKEY_ANNID = "annid";
	public static final String CFGKEY_ANNINFO = "anninfo";
	public static final String CFGKEY_ANNVCFDB = "annvcfdb";
	
	//TSTV
	public static final String CFGKEY_TSTVHOM = "tstvhom";	
	
	//Intervals
	public static final String CFGKEY_INTERX = "interx";
	public static final String CFGKEY_INTERBED = "interbed";
	
	//dbnsfp
	public static final String CFGKEY_DBNSFP = "dbnsfp";
	public static final String CFGKEY_DBNSFPFFIELDS = "dbnsfpFields";
	public static final String CFGKEY_DBNSFPFFIELDSALL = "dbnsfpFieldsall";
	
	
	/**Setting Models**/
	private final SettingsModelString m_snpeff_folder = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_SNPEFF_FOLDER,"");
	private final SettingsModelString m_invcf = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_INVCF,"");
	private final SettingsModelString m_method = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_METHOD,"");
	
	
	/**Filter**/
	private final SettingsModelString m_filterstring = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_FILTERSTRING,"");
	private final SettingsModelDoubleBounded m_filterqual = new SettingsModelDoubleBounded(SnpSiftNodeModel.CFGKEY_FILTERQUAL, 20, 0, Double.MAX_VALUE);
	private final SettingsModelIntegerBounded m_filtercoverage = new SettingsModelIntegerBounded(SnpSiftNodeModel.CFGKEY_FILTERCOVERAGE, 10, 0, Integer.MAX_VALUE);
	private final SettingsModelBoolean m_filterqualbool = new SettingsModelBoolean(CFGKEY_FILTERCOVERAGEBOOL, false);
	private final SettingsModelBoolean m_filtercoveragebool = new SettingsModelBoolean(CFGKEY_FILTERCOVERAGEBOOL, false);

	/**Annotate**/
	private final SettingsModelOptionalString m_anninfo = new SettingsModelOptionalString(
			SnpSiftNodeModel.CFGKEY_ANNINFO,"",false);
	private final SettingsModelBoolean m_annid = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_ANNID, false);
	private final SettingsModelString m_annvcfdb = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_ANNVCFDB,"");
	
	/**TSTV**/
	private final SettingsModelString m_tstvhom = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_TSTVHOM,"");	
	
	/**Intervals**/
	private final SettingsModelString m_interbed = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_INTERBED,"");
	private final SettingsModelBoolean m_interx = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_INTERX, false);

	/**dbnsfp**/
	private final SettingsModelString m_dbnsfp = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_DBNSFP,"");
	private final SettingsModelOptionalString m_dbnsfpfields = new SettingsModelOptionalString(
			SnpSiftNodeModel.CFGKEY_DBNSFPFFIELDS,"",false);
	private final SettingsModelBoolean m_dbnsfpfieldsall = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_DBNSFPFFIELDSALL, false);

	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(SnpSiftNodeModel.class);
	
    protected SnpSiftNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(1, 0);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	/**Initialize logfile**/
    	String logfile = m_snpeff_folder.getStringValue() + "/logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("SnpEff"));
    	/**logfile initialized**/
    	
        String folder = m_snpeff_folder.getStringValue();
        String stdOutFile = "";
        ArrayList<String> command = new ArrayList<String>();
        
        command.add("java");
        command.add("-jar "+folder+"/SnpSift.jar");
        boolean filter = false;
        
        String basename = m_invcf.getStringValue().substring(0,m_invcf.getStringValue().lastIndexOf("."));

        /**FILTER**/
        if(m_method.getStringValue().equals("Filter")){
        	command.add("filter");
        	command.add("--file "+m_invcf.getStringValue());
        	
        	String filterString = "";
        	
        	if(m_filtercoveragebool.getBooleanValue()){
        		filterString+="(DP >= "+m_filtercoverage.getIntValue()+")";
        		filter = true;
        	}
        	if(m_filterqualbool.getBooleanValue()){
        		filterString+=" & (QUAL >= "+m_filterqual.getDoubleValue()+")";
        		filter = true;
        	}
        	if(!(m_filterstring.getStringValue().equals(""))){
        		if(filter){
        			filterString+=" & ";
        		}
        	}
        	filterString+=m_filterstring.getStringValue();
        	command.add(filterString);
        	stdOutFile = basename+"_filtered.vcf";
        }
        
       /**ANNOTATE**/
        if(m_method.getStringValue().equals("Annotate")){
        	command.add("annotate");
        	if(m_annid.getBooleanValue()){
        		command.add("-id");
        	}
        	if(m_anninfo.isEnabled()){
        		command.add(m_anninfo.getStringValue());
        	}
        	command.add(m_annvcfdb.getStringValue());
        	command.add(m_invcf.getStringValue());
        	
        	stdOutFile = basename+"_annotated.vcf";
        }
        
        /**TSTV**/
        if(m_method.getStringValue().equals("TsTv")){
        	command.add("tstv");
        	command.add(m_tstvhom.getStringValue());
        	command.add(m_invcf.getStringValue());
        	
        	stdOutFile = basename+"_tstv.txt";
        }
        
        /**Intervals**/
        if(m_method.getStringValue().equals("Intervals")){
        	command.add("intervals");
        	command.add("-i "+m_invcf.getStringValue());
        	if(m_interx.getBooleanValue()){
        		command.add("-x");
        	}
        	command.add(m_interbed.getStringValue());
        	stdOutFile = basename+"_intervals.vcf";
        }
        
        /**dbnsfp**/
        if(m_method.getStringValue().equals("Annotate with dbnsfp")){
        	command.add("dbnsfp");
        	if(m_dbnsfpfieldsall.getBooleanValue()){
        		command.add("-a");
        	}else{
        		if(m_dbnsfpfields.isActive() && m_dbnsfpfields.isEnabled()){
        			command.add("-f "+m_dbnsfpfields.getStringValue());
        		}	
        	}
        	command.add(m_dbnsfp.getStringValue());
        	command.add(m_invcf.getStringValue());
        	stdOutFile = basename+"_dbnsfp.vcf";
        }
        

    	/**
    	 * Execute
    	 */
    	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER,stdOutFile);
    	logBuffer.append(ShowOutput.getNodeEndTime());
    	ShowOutput.writeLogFile(logBuffer);
         	
        return null;
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
    	m_filtercoverage.setEnabled(false);
    	m_filterqual.setEnabled(false);
    	
        // TODO: generated method stub
        return null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	m_invcf.saveSettingsTo(settings);
    	m_method.saveSettingsTo(settings);
    	m_snpeff_folder.saveSettingsTo(settings);
         m_filtercoverage.saveSettingsTo(settings);
         m_filterqual.saveSettingsTo(settings);
         m_filterstring.saveSettingsTo(settings);
         m_filtercoveragebool.saveSettingsTo(settings);
         m_filterqualbool.saveSettingsTo(settings);
         m_annid.saveSettingsTo(settings);
         m_anninfo.saveSettingsTo(settings);
         m_annvcfdb.saveSettingsTo(settings);
         m_tstvhom.saveSettingsTo(settings);
         m_interbed.saveSettingsTo(settings);
         m_interx.saveSettingsTo(settings);
         m_dbnsfp.saveSettingsTo(settings);
         m_dbnsfpfields.saveSettingsTo(settings);
         m_dbnsfpfieldsall.saveSettingsTo(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_filtercoveragebool.loadSettingsFrom(settings);
        m_filterqualbool.loadSettingsFrom(settings);
        m_filtercoverage.loadSettingsFrom(settings);
        m_filterqual.loadSettingsFrom(settings);
        m_filterstring.loadSettingsFrom(settings);
        m_invcf.loadSettingsFrom(settings);
        m_method.loadSettingsFrom(settings);
        m_snpeff_folder.loadSettingsFrom(settings);
        m_annid.loadSettingsFrom(settings);
        m_anninfo.loadSettingsFrom(settings);
        m_annvcfdb.loadSettingsFrom(settings);
        m_tstvhom.loadSettingsFrom(settings);
        m_interbed.loadSettingsFrom(settings);
        m_interx.loadSettingsFrom(settings);
        m_dbnsfp.loadSettingsFrom(settings);
        m_dbnsfpfields.loadSettingsFrom(settings);
        m_dbnsfpfieldsall.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        m_filtercoverage.validateSettings(settings);
        m_filterqual.validateSettings(settings);
        m_filtercoveragebool.validateSettings(settings);
        m_filterqualbool.validateSettings(settings);
        m_filterstring.validateSettings(settings);
        m_invcf.validateSettings(settings);
        m_method.validateSettings(settings);
        m_snpeff_folder.validateSettings(settings);
        m_annid.validateSettings(settings);
        m_anninfo.validateSettings(settings);
        m_annvcfdb.validateSettings(settings);
        m_tstvhom.validateSettings(settings);
        m_interx.validateSettings(settings);
        m_interbed.validateSettings(settings);
        m_dbnsfp.validateSettings(settings);
        m_dbnsfpfields.validateSettings(settings);
        m_dbnsfpfieldsall.validateSettings(settings);
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


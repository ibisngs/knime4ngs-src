package de.helmholtz_muenchen.ibis.ngs.snpsift;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;



/**
 * This is the model implementation of SnpSift.
 * 
 *
 * @author Max
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
        String com = "java -jar "+folder+"/SnpSift.jar "; 
        boolean filter = false;
        
        String basename = m_invcf.getStringValue().substring(0,m_invcf.getStringValue().lastIndexOf("."));

        /**FILTER**/
        if(m_method.getStringValue().equals("Filter")){
        	com+="filter --file "+m_invcf.getStringValue()+" \"";
        	if(m_filtercoveragebool.getBooleanValue()){
        		com+=" (DP >= "+m_filtercoverage.getIntValue()+")";
        		filter = true;
        	}
        	if(m_filterqualbool.getBooleanValue()){
        		com+=" & (QUAL >= "+m_filterqual.getDoubleValue()+")";
        		filter = true;
        	}
        	if(!(m_filterstring.getStringValue().equals(""))){
        		if(filter){
        			com+=" & ";
        		}
        	}
        	com+=m_filterstring.getStringValue()+"\"";
        	com+=" > "+basename+"_filtered.vcf";
        }
        
       /**ANNOTATE**/
        if(m_method.getStringValue().equals("Annotate")){
        	com+="annotate ";
        	if(m_annid.getBooleanValue()){
        		com+="-id ";
        	}
        	if(m_anninfo.isEnabled()){
        		com+=m_anninfo.getStringValue()+" ";
        	}
        	com+=m_annvcfdb.getStringValue()+" "+m_invcf.getStringValue()+" > "+basename+"_annotated.vcf";
        }
        
        /**TSTV**/
        if(m_method.getStringValue().equals("TsTv")){
        	com+="tstv "+m_tstvhom.getStringValue()+" "+m_invcf.getStringValue()+" > "+basename+"_tstv.txt";
        }
        
        /**Intervals**/
        if(m_method.getStringValue().equals("Intervals")){
        	com+="intervals -i "+m_invcf.getStringValue()+" ";
        	if(m_interx.getBooleanValue()){
        		com+="-x ";
        	}
        	com+=m_interbed.getStringValue() + " > "+basename+"_intervals.vcf";
        }
        
        /**dbnsfp**/
        if(m_method.getStringValue().equals("Annotate with dbnsfp")){
        	com+="dbnsfp ";
        	if(m_dbnsfpfieldsall.getBooleanValue()){
        		com+="-a ";
        	}else{
        		if(m_dbnsfpfields.isActive() && m_dbnsfpfields.isEnabled()){
        			com+="-f "+m_dbnsfpfields.getStringValue()+" ";
        		}	
        	}
        	com+=m_dbnsfp.getStringValue()+" "+m_invcf.getStringValue()+" > "+basename+"_dbnsfp.vcf";
        }
        
    	System.out.println(com);
    	System.out.println(logfile);
    	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", com);
    	Process p1 = b.start();
    	p1.waitFor();
        logBuffer.append(ShowOutput.getLogEntry(p1, com));
        ShowOutput.writeLogFile(new StringBuffer(ShowOutput.getNodeEndTime()));
         	
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


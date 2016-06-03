package de.helmholtz_muenchen.ibis.ngs.gatkvariantfiltration;


import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;


/**
 * This is the model implementation of GATKVariantFiltration.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKVariantFiltrationNodeModel extends GATKNodeModel {
    
	
	/**
	 * Config Keys
	 */
//	public static final String CFGKEY_QUAL = "Quality score";
	public static final String CFGKEY_QD = "QualByDepth";
	public static final String CFGKEY_FS = "Fisherstrand";
	public static final String CFGKEY_MQ = "MappingQuality";
	public static final String CFGKEY_HS = "HaplotypeScore";
	public static final String CFGKEY_MQR = "MappingQualityRankSum";
	public static final String CFGKEY_RPR = "ReadPosRankSum";
	public static final String CFGKEY_INFOFilterString = "INFOFilterString";
	public static final String CFGKEY_INFOFilterName = "INFOFilterName";
	
	public static final String CFGKEY_DP = "DepthOfReads";
	public static final String CFGKEY_GQ = "GenotypeQuality";
	public static final String CFGKEY_NOCALL = "NoCall";
	public static final String CFGKEY_FORMATFilterString = "FORMATFilterString";
	public static final String CFGKEY_FORMATFilterName = "FORMATFilterName";
	
	/**
	 * The SettingsModels
	 */
//	private final SettingsModelOptionalString m_QUAL= new SettingsModelOptionalString(CFGKEY_QUAL, "<50.0",true);
	private final SettingsModelOptionalString m_QD= new SettingsModelOptionalString(CFGKEY_QD, "<2.0",false);
	private final SettingsModelOptionalString m_FS= new SettingsModelOptionalString(CFGKEY_FS, ">60.0",false);
	private final SettingsModelOptionalString m_MQ= new SettingsModelOptionalString(CFGKEY_MQ, "<40.0",false);
	private final SettingsModelOptionalString m_HS= new SettingsModelOptionalString(CFGKEY_HS, ">13.0",false);
	private final SettingsModelOptionalString m_MQR= new SettingsModelOptionalString(CFGKEY_MQR, "<-12.5",false);
	private final SettingsModelOptionalString m_RPR= new SettingsModelOptionalString(CFGKEY_RPR, "<-8.0",false);
	private final SettingsModelOptionalString m_INFOFilterString = new SettingsModelOptionalString(CFGKEY_INFOFilterString, "Value1<X||Value2>Y||...",false);
	private final SettingsModelOptionalString m_INFOFilterName = new SettingsModelOptionalString(CFGKEY_INFOFilterName, "GATKVariantFiltration",false);
	
	private final SettingsModelOptionalString m_DP= new SettingsModelOptionalString(CFGKEY_DP, "<8.0",true);
	private final SettingsModelOptionalString m_GQ= new SettingsModelOptionalString(CFGKEY_GQ, "<20.0",true);
	private final SettingsModelBoolean m_NOCALL = new SettingsModelBoolean(CFGKEY_NOCALL, true);
	private final SettingsModelOptionalString m_FORMATFilterString = new SettingsModelOptionalString(CFGKEY_FORMATFilterString, "Value1<X||Value2>Y||...",false);
	private final SettingsModelOptionalString m_FORMATFilterName = new SettingsModelOptionalString(CFGKEY_FORMATFilterName, "GATKVariantFiltration",false);
	
	//The Output Col Names
	public static final String OUT_COL1 = "FILTERED_VARIANTS";
	
	private String OUTFILE;
	private int vcf_index;
	
    /**
     * Constructor for the node model.
     */
    protected GATKVariantFiltrationNodeModel() {
        super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
        
        addSetting(m_DP);
        addSetting(m_GQ);
        addSetting(m_FS);
        addSetting(m_HS);
        addSetting(m_MQ);
        addSetting(m_MQR);
        addSetting(m_QD);
//        addSetting(m_QUAL);
        addSetting(m_RPR);
        addSetting(m_INFOFilterString);
        addSetting(m_INFOFilterName);
        addSetting(m_FORMATFilterString);
        addSetting(m_FORMATFilterName);
        addSetting(m_NOCALL);
    }

    /**
     * Adds the values of SettingsModelString to the filterString if model is enabled
     * @param toAdd
     */
    private StringBuffer addToFilterString(SettingsModelOptionalString toAdd, String FieldName, StringBuffer filterString){
    	if(toAdd.isActive()){
        	if(filterString.length()==0 && !(FieldName.equals(""))){
        		filterString.append(FieldName+toAdd.getStringValue());
        	}else if(filterString.length() > 0 && !(FieldName.equals(""))){
        		filterString.append("||"+FieldName+toAdd.getStringValue());
        	}else{
        		filterString.append("||"+toAdd.getStringValue());
        	}
    	}
    	return filterString;
    }

	@Override
	protected String getCommandParameters(BufferedDataTable[] inData)
			throws InvalidSettingsException {
		/**
    	 * Check INFILE
    	 */
    	String INFILE;
    	try{
    		INFILE = inData[0].iterator().next().getCell(vcf_index).toString();
    		if(!INFILE.endsWith(".vcf")){
    			throw new InvalidSettingsException("A cell of the input table has to be the path to VCF infile but it is "+INFILE);
    		}
    	}catch(IndexOutOfBoundsException e){
    			throw new InvalidSettingsException("A cell of the input table has to be the path to VCF infile but it is empty.");
    	}
    	
    	ArrayList<String> command = new ArrayList<String>();

    	
    	command.add("-V "+INFILE);
    	
    	OUTFILE = IO.replaceFileExtension(INFILE, ".VariantFiltration.vcf");
    	
    	/**
    	 * String that holds the complete filter options
    	 */
    	StringBuffer filterStringINFOFIELD = new StringBuffer();
    	StringBuffer filterStringFORMATFIELD = new StringBuffer();
    	
    	/**
    	 * Create Filter String for INFO Field
    	 */
//    	filterStringINFOFIELD = addToFilterString(m_QUAL,"QUAL",filterStringINFOFIELD);
    	filterStringINFOFIELD = addToFilterString(m_QD,"QD",filterStringINFOFIELD);
    	filterStringINFOFIELD = addToFilterString(m_FS,"FS",filterStringINFOFIELD);
    	filterStringINFOFIELD = addToFilterString(m_MQ,"MQ",filterStringINFOFIELD);
    	filterStringINFOFIELD = addToFilterString(m_HS,"HaplotypeScore",filterStringINFOFIELD);
    	filterStringINFOFIELD = addToFilterString(m_MQR,"MappingQualityRankSum",filterStringINFOFIELD);
    	filterStringINFOFIELD = addToFilterString(m_RPR,"ReadPosRankSum",filterStringINFOFIELD);
    	filterStringINFOFIELD = addToFilterString(m_INFOFilterString,"",filterStringINFOFIELD);
		
    	if(filterStringINFOFIELD.length()!=0){
			command.add("--filterExpression");
			command.add(filterStringINFOFIELD.toString());
		}
    	
		/**
		 * Create Filter String for FORMAT Field
		 */
		filterStringFORMATFIELD = addToFilterString(m_DP,"DP",filterStringFORMATFIELD);
		filterStringFORMATFIELD = addToFilterString(m_GQ,"GQ",filterStringFORMATFIELD);
    	filterStringFORMATFIELD = addToFilterString(m_FORMATFilterString,"",filterStringFORMATFIELD);

		if(filterStringFORMATFIELD.length()!=0){
			command.add("--genotypeFilterExpression");    	
    		command.add(filterStringFORMATFIELD.toString());
    	}
    	
    	
    	if(filterStringINFOFIELD.length()!=0){
        	command.add("--filterName");
        	if(m_INFOFilterName.isActive()){
        		command.add(m_INFOFilterName.getStringValue());
        	} else {
        		command.add("GATKVariantFiltration");
        	}
    	}
    	
    	if(filterStringFORMATFIELD.length()!=0){
    		command.add("--genotypeFilterName");
    		if(m_FORMATFilterName.isActive()) {
        		command.add(m_FORMATFilterName.getStringValue());
    		} else {
    			command.add("GATKVariantFiltration");
    		}
    		
    		if(m_NOCALL.getBooleanValue()) {
        		command.add("--setFilteredGtToNocall");
        	}
    	}
    	
    	return StringUtils.join(command, " ");
	}

	@Override
	protected String getCommandWalker() {
		return "VariantFiltration";
	}


	@Override
	protected String getOutfile() {
		return OUTFILE;
	}

	@Override
	protected boolean checkInputCellType(DataTableSpec[] inSpecs) {
		vcf_index = CompatibilityChecker.getFirstIndexCellType(inSpecs[0], "VCFCell");
		return (vcf_index>-1);
	}

	@Override
	protected DataType getOutColType() {
		return VCFCell.TYPE;
	}

	@Override
	protected void extraConfig() throws InvalidSettingsException {
		// TODO Auto-generated method stub
		
	}
}


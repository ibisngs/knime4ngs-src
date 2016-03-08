package de.helmholtz_muenchen.ibis.ngs.gatkMergeTwoVCFs;

import java.io.File;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.port.PortType;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of MergeTwoVCFs.
 * 
 *
 * @author Kaarin Ahomaa
 */
public class MergeTwoVCFsNodeModel extends GATKNodeModel {
    
	protected MergeTwoVCFsNodeModel(PortType[] INPORTS, PortType[] OUTPORTS) {
		super(INPORTS, OUTPORTS);
		
	}
	
	//configuration keys
	static final String CFGKEY_INPUT1_TAG = "input1_tag";
	static final String CFGKEY_INPUT2_TAG = "input2_tag";
	static final String CFGKEY_PRIORITIZE = "priority";
	static final String CFGKEY_OUTFOLDER = "outfolder";
	static final String CFGKEY_FILTEREDRECORDSMERGETYPE = "filteredrecordsmergetype";
	

	
	//settings models
	public static final String CFGKEY_GENOTYPEMERGEOPTION 		= "genotypemergeoption";
	
	private final SettingsModelString m_GENOTYPEMERGEOPTION	= new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_GENOTYPEMERGEOPTION, "");
	private final SettingsModelString m_INPUT1_TAG = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_INPUT1_TAG, "");
	private final SettingsModelString m_INPUT2_TAG = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_INPUT2_TAG,"");
	private final SettingsModelString m_PRIORITIZE = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_PRIORITIZE,"");
	private final SettingsModelString m_OUTFOLDER = new SettingsModelString(CFGKEY_OUTFOLDER, "");
	private final SettingsModelString m_FILTEREDRECORDSMERGETYPE = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_FILTEREDRECORDSMERGETYPE, "");
	
	
	private String outfile;
	private int vcf_ind1, vcf_ind2;

    /**
     * Constructor for the node model.
     */
    protected MergeTwoVCFsNodeModel() {

    	 //Specify the amount of input and output ports needed.
    	super(OptionalPorts.createOPOs(2), OptionalPorts.createOPOs(1));
    	addSetting(m_GENOTYPEMERGEOPTION);
		addSetting(m_INPUT1_TAG);
		addSetting(m_INPUT2_TAG);
		addSetting(m_PRIORITIZE);
		addSetting(m_OUTFOLDER);
    }

    
	@Override
	protected String getCommandParameters(BufferedDataTable[] inData) throws InvalidSettingsException {
		
		//input file

		String  vcf1 = inData[0].iterator().next().getCell(vcf_ind1).toString();
		String  vcf2 = inData[1].iterator().next().getCell(vcf_ind2).toString();
		ArrayList<String> command 	= new ArrayList<String>();
		
		outfile = m_OUTFOLDER.getStringValue()+ System.getProperty("file.separator")+ new File(vcf1).getName();
		outfile = IO.replaceFileExtension(outfile, "MERGED.vcf"); 
		
		String filterOption = m_FILTEREDRECORDSMERGETYPE.getStringValue();
		command.add("--filteredrecordsmergetype "+filterOption);
		
		
		String mergeOption = m_GENOTYPEMERGEOPTION.getStringValue();
		command.add("--genotypemergeoption "+mergeOption);
		if(mergeOption.equals("PRIORITIZE")) {
			command.add("--variant:"+m_INPUT1_TAG.getStringValue()+" "+vcf1);
			command.add("--variant:"+m_INPUT2_TAG.getStringValue()+ " "+vcf2);
			command.add("-priority "+m_PRIORITIZE.getStringValue());
		} else {
			command.add("--variant "+vcf1);
			command.add("--variant "+vcf2);
		}
		return StringUtils.join(command, " ");
	}
		

	@Override
	protected String getCommandWalker() {
		return "CombineVariants";
	}
	

	
	@Override
	protected boolean checkInputCellType(DataTableSpec[] inSpecs) {
		vcf_ind1 = -1;
		
		for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
			if(inSpecs[0].getColumnSpec(i).getType().toString().equals("VCFCell")) {
				vcf_ind1 = i;
			}
		}
		
		vcf_ind2 = -1;
		
		for(int i = 0; i < inSpecs[1].getNumColumns(); i++) {
			if(inSpecs[1].getColumnSpec(i).getType().toString().equals("VCFCell")) {
				vcf_ind2 = i;
			}
		} 
		
		
		return (vcf_ind1>-1 && vcf_ind2>-1);
	}

	
	@Override
	protected DataType getOutColType() {
		return VCFCell.TYPE;
	}


	@Override
	protected void extraConfig() throws InvalidSettingsException {
		String outfolder_warning = CheckUtils.checkDestinationDirectory(m_OUTFOLDER.getStringValue());
		if(outfolder_warning!=null) {
			setWarningMessage(outfolder_warning);
		}
	}


//	@Override
//	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
//		m_GENOTYPEMERGEOPTION.saveSettingsTo(settings);
//		m_INPUT1_TAG.saveSettingsTo(settings);
//		m_INPUT2_TAG.saveSettingsTo(settings);
//		m_PRIORITIZE.saveSettingsTo(settings);
//		m_OUTFOLDER.saveSettingsTo(settings);
//		
//	}
//
//
//	@Override
//	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
//			throws InvalidSettingsException {
//		m_GENOTYPEMERGEOPTION.loadSettingsFrom(settings);
//		m_INPUT1_TAG.loadSettingsFrom(settings);
//		m_INPUT2_TAG.loadSettingsFrom(settings);
//		m_PRIORITIZE.loadSettingsFrom(settings);
//		m_OUTFOLDER.loadSettingsFrom(settings);
//	}
//
//
//	@Override
//	protected void validateExtraSettings(NodeSettingsRO settings)
//			throws InvalidSettingsException {
//		m_GENOTYPEMERGEOPTION.validateSettings(settings);
//		m_INPUT1_TAG.validateSettings(settings);
//		m_INPUT2_TAG.validateSettings(settings);
//		m_PRIORITIZE.validateSettings(settings);
//		m_OUTFOLDER.validateSettings(settings);
//	}


	@Override
	protected String getOutfile() {
		return outfile;
	}

}

package de.helmholtz_muenchen.ibis.ngs.gatkMergeTwoVCFs;

import java.io.File;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.port.PortType;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
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
	
	//settings models
	public static final String CFGKEY_GENOTYPEMERGEOPTION 		= "genotypemergeoption";
	
	private final SettingsModelString m_GENOTYPEMERGEOPTION	= new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_GENOTYPEMERGEOPTION, "");
	private final SettingsModelString m_INPUT1_TAG = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_INPUT1_TAG, "");
	private final SettingsModelString m_INPUT2_TAG = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_INPUT2_TAG,"");
	private final SettingsModelString m_PRIORITIZE = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_PRIORITIZE,"");
	
	private String OUTFILE, LOCKFILE;
	private int vcf_ind1, vcf_ind2;

    /**
     * Constructor for the node model.
     */
    protected MergeTwoVCFsNodeModel() {

    	 //Specify the amount of input and output ports needed.
    	super(OptionalPorts.createOPOs(2), OptionalPorts.createOPOs(1));
    }

    
	@Override
	protected String getCommandParameters(BufferedDataTable[] inData) throws InvalidSettingsException {
		
		//input file

		String  vcf1 = inData[0].iterator().next().getCell(vcf_ind1).toString();
		String  vcf2 = inData[1].iterator().next().getCell(vcf_ind2).toString();
		ArrayList<String> command 	= new ArrayList<String>();
		
		this.OUTFILE = IO.replaceFileExtension(vcf1, ".MERGEDVCF.vcf");
		this.LOCKFILE = IO.replaceFileExtension(vcf1, SuccessfulRunChecker.LOCK_ENDING);
			

		
		
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
	protected File getLockFile() {
		return new File(this.LOCKFILE);
	}
	
	
	@Override
	protected String getOutfile() {
		return this.OUTFILE;
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
		
		for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
			if(inSpecs[0].getColumnSpec(i).getType().toString().equals("VCFCell")) {
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
	protected void extraConfig()
			throws InvalidSettingsException {
		
	}


	@Override
	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
		m_GENOTYPEMERGEOPTION.saveSettingsTo(settings);
		m_INPUT1_TAG.saveSettingsTo(settings);
		m_INPUT2_TAG.saveSettingsTo(settings);
		m_PRIORITIZE.saveSettingsTo(settings);
		
	}


	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_GENOTYPEMERGEOPTION.loadSettingsFrom(settings);
		m_INPUT1_TAG.loadSettingsFrom(settings);
		m_INPUT2_TAG.loadSettingsFrom(settings);
		m_PRIORITIZE.loadSettingsFrom(settings);
	}


	@Override
	protected void validateExtraSettings(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_GENOTYPEMERGEOPTION.validateSettings(settings);
		m_INPUT1_TAG.validateSettings(settings);
		m_INPUT2_TAG.validateSettings(settings);
		m_PRIORITIZE.validateSettings(settings);
	}

}

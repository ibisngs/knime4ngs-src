package de.helmholtz_muenchen.ibis.ngs.gatkMergeTwoVCFs;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataRow;
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
	static final String CFGKEY_INPUT1 = "input1.vcf";
	static final String CFGKEY_INPUT2 = "input2.vcf";
	static final String CFGKEY_PRIORITIZE = "priority";
	
	//settings models
	public static final String CFGKEY_GENOTYPEMERGEOPTION 		= "genotypemergeoption";
	
	private final SettingsModelString m_GENOTYPEMERGEOPTION	= new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_GENOTYPEMERGEOPTION, "");
	private final SettingsModelString m_INPUT1 = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_INPUT1, "");
	private final SettingsModelString m_INPUT2 = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_INPUT2,"");
	private final SettingsModelString m_PRIORITIZE = new SettingsModelString(MergeTwoVCFsNodeModel.CFGKEY_PRIORITIZE,"");
	
	private String OUTFILE, LOCKFILE;
	private int vcf_index;
	
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
		String vcf_infile1 = inData[0].iterator().next().getCell(vcf_index).toString();
		String vcf_infile2 = inData[0].iterator().next().getCell(vcf_index).toString();
    	
    	
		Iterator <DataRow> it = inData[0].iterator();
		ArrayList<String> command 	= new ArrayList<String>();
		boolean first = true;
		while(it.hasNext()){
			DataRow row = it.next();
			String INFILE = row.getCell(vcf_index).toString();
			
			if(first){
				this.OUTFILE = IO.replaceFileExtension(INFILE, ".MERGEDVCF.vcf");
				this.LOCKFILE = IO.replaceFileExtension(INFILE, SuccessfulRunChecker.LOCK_ENDING);
				first=false;
			}

			command.add("--variant "+m_INPUT1.getStringValue());
			command.add("--variant "+m_INPUT2.getStringValue());
		}
		
		command.add("--genotypemergeoption "+m_GENOTYPEMERGEOPTION.getStringValue());
		command.add("-priority "+m_PRIORITIZE.getStringValue());
	
		
		
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
		vcf_index = -1;
		
		for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
			if(inSpecs[0].getColumnSpec(i).getType().toString().equals("VCFCell")) {
				vcf_index = i;
			}
		}
		return (vcf_index>-1);
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
		m_INPUT1.saveSettingsTo(settings);
		m_INPUT2.saveSettingsTo(settings);
		m_PRIORITIZE.saveSettingsTo(settings);
		
	}


	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_GENOTYPEMERGEOPTION.loadSettingsFrom(settings);
		m_INPUT1.loadSettingsFrom(settings);
		m_INPUT2.loadSettingsFrom(settings);
		m_PRIORITIZE.loadSettingsFrom(settings);
	}


	@Override
	protected void validateExtraSettings(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_GENOTYPEMERGEOPTION.validateSettings(settings);
		m_INPUT1.validateSettings(settings);
		m_INPUT2.validateSettings(settings);
		m_PRIORITIZE.validateSettings(settings);
	}

}

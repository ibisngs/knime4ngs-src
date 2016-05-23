package de.helmholtz_muenchen.ibis.ngs.gatkCombineVCFs;


import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.ngs.vcfmerger.VCFMergerNodeModel;
import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.IO;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of CombineVCFs.
 * 
 *
 * @author Kaarin Ahomaa
 */
public class CombineVCFsNodeModel extends GATKNodeModel {
    
	public static final String CFGKEY_GENOTYPEMERGEOPTION = "genotypemergeoption";
	static final String CFGKEY_OUTFOLDER = "outfolder";
	
	private final SettingsModelString m_GENOTYPEMERGEOPTION	= new SettingsModelString(VCFMergerNodeModel.CFGKEY_GENOTYPEMERGEOPTION, "");
	private final SettingsModelString m_OUTFOLDER = new SettingsModelString(VCFMergerNodeModel.CFGKEY_OUTFOLDER, "");
	
	private String outfile;
	private int vcf_index;
	
    /**
     * Constructor for the node model.
     */
    protected CombineVCFsNodeModel() {
    
    	super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
    	addSetting(m_GENOTYPEMERGEOPTION);
		addSetting(m_OUTFOLDER);
    }

    
	@Override
	protected String getCommandParameters(BufferedDataTable[] inData) throws InvalidSettingsException {
		
		
		
		
		/**
    	 * Check INFILE
    	 */
		
		Iterator <DataRow> it = inData[0].iterator();
		ArrayList<String> command 	= new ArrayList<String>();
		boolean first = true;
		while(it.hasNext()){
			DataRow row = it.next();
			String INFILE = row.getCell(vcf_index).toString();
			
			
			if(first){
				outfile = m_OUTFOLDER.getStringValue()+ System.getProperty("file.separator")+ new File(INFILE).getName();
				outfile = IO.replaceFileExtension(outfile, ".ALLVARIANTS.vcf");
				first=false;
			}

			command.add("--variant "+INFILE);
		}
		
		command.add("--genotypemergeoption "+m_GENOTYPEMERGEOPTION.getStringValue());
		
		return StringUtils.join(command, " ");
	}
		
	

	@Override
	protected String getCommandWalker() {
		return "CombineVariants";
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
	protected void extraConfig() throws InvalidSettingsException {
		String outfolder_warning = CheckUtils.checkDestinationDirectory(m_OUTFOLDER.getStringValue());
		if(outfolder_warning!=null) {
			setWarningMessage(outfolder_warning);
		}
	}

	

	 
	
	
//	@Override
//	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
//		m_GENOTYPEMERGEOPTION.saveSettingsTo(settings);
//		m_OUTFOLDER.saveSettingsTo(settings);
//		
//	}
//
//
//	@Override
//	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
//			throws InvalidSettingsException {
//		m_GENOTYPEMERGEOPTION.loadSettingsFrom(settings);
//		m_OUTFOLDER.loadSettingsFrom(settings);
//	}
//
//
//	@Override
//	protected void validateExtraSettings(NodeSettingsRO settings)
//			throws InvalidSettingsException {
//		m_GENOTYPEMERGEOPTION.validateSettings(settings);
//		m_OUTFOLDER.validateSettings(settings);
//	}

	@Override
	protected String getOutfile() {
		return outfile;
	}

}


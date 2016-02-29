package de.helmholtz_muenchen.ibis.ngs.gatkgenotypegvcfs;



import java.io.File;
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
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GenotypeGVCFs.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GenotypeGVCFsNodeModel extends GATKNodeModel {
    
	static final String CFGKEY_NT_FILE = "NT";
	static final String CFGKEY_OUTFOLDER = "outfolder";

	private final SettingsModelIntegerBounded m_NT = new SettingsModelIntegerBounded(GenotypeGVCFsNodeModel.CFGKEY_NT_FILE, 1, 1, Integer.MAX_VALUE);
	private final SettingsModelString m_OUTFOLDER = new SettingsModelString(CFGKEY_OUTFOLDER, "");
	
	private String OUTFILE;
	private int gvcf_index;

    protected GenotypeGVCFsNodeModel(int INPORTS, int OUTPORTS) {
		super(OptionalPorts.createOPOs(INPORTS), OptionalPorts.createOPOs(OUTPORTS));
   	}

	@Override
	protected String getCommandParameters(BufferedDataTable[] inData) throws InvalidSettingsException {
		
		Iterator <DataRow> it = inData[0].iterator();
		ArrayList<String> command 	= new ArrayList<String>();
		boolean first = true;
		while(it.hasNext()){
			DataRow row = it.next();
			String INFILE = row.getCell(gvcf_index).toString();
			
			if(first){
				OUTFILE = m_OUTFOLDER.getStringValue()+ System.getProperty("file.separator")+ new File(INFILE).getName();
				OUTFILE = IO.replaceFileExtension(OUTFILE, ".GenotypedVariants.vcf");
				first=false;
			}

			command.add("--variant "+INFILE);

		}

		command.add("-nt "+m_NT.getIntValue());
		return StringUtils.join(command, " ");
	}

	@Override
	protected String getCommandWalker() {
		return "GenotypeGVCFs";
	}


	@Override
	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
		m_NT.saveSettingsTo(settings);
		m_OUTFOLDER.saveSettingsTo(settings);
		
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_NT.loadSettingsFrom(settings);
		m_OUTFOLDER.loadSettingsFrom(settings);
		
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_NT.validateSettings(settings);
		m_OUTFOLDER.validateSettings(settings);
	}
	
	@Override
	protected boolean checkInputCellType(DataTableSpec[] inSpecs) {
		gvcf_index = -1;
		
		for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
    		if(inSpecs[0].getColumnSpec(i).getType().toString().equals("GVCFCell")) {
    			gvcf_index = i;
    		}
    	}
		return (gvcf_index>-1);
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
	
	@Override
	protected String getOutfile() {
		return OUTFILE;
	}
}


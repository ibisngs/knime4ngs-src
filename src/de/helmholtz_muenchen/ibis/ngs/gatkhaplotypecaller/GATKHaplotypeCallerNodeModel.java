package de.helmholtz_muenchen.ibis.ngs.gatkhaplotypecaller;


import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;

/**
 * This is the model implementation of GATKHaplotypeCaller.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKHaplotypeCallerNodeModel extends GATKNodeModel {

	
	public static final String CFGKEY_BED_FILE 			= "BED_FILE";
	public static final String CFGKEY_BED_FILE_CHECKBOX = "BED_FILE_CHECKBOX";
	
	private final SettingsModelString m_BED_FILE = new SettingsModelString(GATKHaplotypeCallerNodeModel.CFGKEY_BED_FILE, "");
    private final SettingsModelBoolean m_BED_FILE_CHECKBOX = new SettingsModelBoolean(GATKHaplotypeCallerNodeModel.CFGKEY_BED_FILE_CHECKBOX, false);
	private String OUTFILE; 
	
	protected GATKHaplotypeCallerNodeModel(int INPORTS, int OUTPORTS) {
		super(INPORTS, OUTPORTS);
	}

	@Override
	protected String getCommandParameters(final BufferedDataTable[] inData) {
		
		String INFILE				= inData[0].iterator().next().getCell(0).toString();
		ArrayList<String> command 	= new ArrayList<String>();
		command.add("-I "+INFILE);
		command.add("--emitRefConfidence GVCF");
		command.add("--variant_index_type LINEAR");
		command.add("--variant_index_parameter 128000");
		command.add("-L "+m_BED_FILE.getStringValue());
		
		this.OUTFILE = IO.replaceFileExtension(INFILE, ".gvcf");

		return StringUtils.join(command, " ");
	}

	@Override
	protected String getCommandWalker() {
		return "HaplotypeCaller";
	}

	@Override
	protected String getOutfile() {
		return this.OUTFILE;
	}

	@Override
	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
		m_BED_FILE.saveSettingsTo(settings);
		m_BED_FILE_CHECKBOX.saveSettingsTo(settings);
		
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_BED_FILE.loadSettingsFrom(settings);
		m_BED_FILE_CHECKBOX.loadSettingsFrom(settings);
		
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_BED_FILE.validateSettings(settings);
		m_BED_FILE_CHECKBOX.validateSettings(settings);
		
	}
    


}


package de.helmholtz_muenchen.ibis.ngs.gatkhaplotypecaller;


import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataRow;
import org.knime.core.data.container.CloseableRowIterator;
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
	
//	private ArrayList<String> INFILE_arrList = new ArrayList<>(3);
	private final SettingsModelString m_BED_FILE = new SettingsModelString(GATKHaplotypeCallerNodeModel.CFGKEY_BED_FILE, "");
    private final SettingsModelBoolean m_BED_FILE_CHECKBOX = new SettingsModelBoolean(GATKHaplotypeCallerNodeModel.CFGKEY_BED_FILE_CHECKBOX, false);
	private String INFILE; 
    private String OUTFILE; 
	
	protected GATKHaplotypeCallerNodeModel(int INPORTS, int OUTPORTS) {
		super(INPORTS, OUTPORTS);
	}

	@Override
	protected String getCommandParameters(final BufferedDataTable[] inData) {
		
		CloseableRowIterator it = inData[0].iterator();
		while (it.hasNext()) {
			DataRow row = it.next();
//			INFILE_arrList.add(row.getCell(0).toString());
			INFILE = row.getCell(0).toString();
		}
//		String INFILE				= inData[0].iterator().next().getCell(0).toString();
		ArrayList<String> command 	= new ArrayList<String>();
//		for(String s : INFILE_arrList) {
//    		command.add("-I "+s);
//
//		}
		command.add("-I "+INFILE);
//		if (INFILE_arrList.size()==1) {
		command.add("--emitRefConfidence GVCF");
//		}
		command.add("--variant_index_type LINEAR");
		command.add("--variant_index_parameter 128000");
		if(m_BED_FILE_CHECKBOX.isEnabled())
			command.add("-L "+m_BED_FILE.getStringValue());
		
		this.OUTFILE = IO.replaceFileExtension(INFILE/*_arrList.get(0)*/, ".gvcf");

		String commandLine = StringUtils.join(command, " ");
		
		String reffile = "";
		String[] lines = commandLine.split(" ");
		for (int i = 0; i < lines.length-1; i++) {
			if (lines[i] == "-R") {
				reffile = lines[i+1];
//				System.out.println(reffile);
				break;
			}
		}
		
		/**
    	 * push the reference file extra to flow variable for the phasers in next step
    	 */
    	pushFlowVariableString("Reference", reffile); 
    	/**
    	 * 
    	 */
		return commandLine;
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


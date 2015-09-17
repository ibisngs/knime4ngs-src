package de.helmholtz_muenchen.ibis.ngs.gatkhaplotypecaller;


import java.io.File;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataRow;
import org.knime.core.data.container.CloseableRowIterator;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GATKHaplotypeCaller.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKHaplotypeCallerNodeModel extends GATKNodeModel {

	private String INFILE; 
    private String OUTFILE;
    private String LOCKFILE;
	
	protected GATKHaplotypeCallerNodeModel(int INPORTS, int OUTPORTS) {
		super(OptionalPorts.createOPOs(INPORTS), OptionalPorts.createOPOs(OUTPORTS));
	}

	@Override
	protected String getCommandParameters(final BufferedDataTable[] inData) {
		
		CloseableRowIterator it = inData[0].iterator();
		while (it.hasNext()) {
			DataRow row = it.next();
			INFILE = row.getCell(0).toString();
		}
		ArrayList<String> command 	= new ArrayList<String>();

		command.add("-I "+INFILE);
		command.add("--emitRefConfidence GVCF");
		command.add("--variant_index_type LINEAR");
		command.add("--variant_index_parameter 128000");
		
		this.OUTFILE = IO.replaceFileExtension(INFILE/*_arrList.get(0)*/, ".gvcf");
		this.LOCKFILE = IO.replaceFileExtension(INFILE, SuccessfulRunChecker.LOCK_ENDING);
		String commandLine = StringUtils.join(command, " ");
		
		String reffile = "";
		String[] lines = commandLine.split(" ");
		for (int i = 0; i < lines.length-1; i++) {
			if (lines[i] == "-R") {
				reffile = lines[i+1];
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
		
		
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
			throws InvalidSettingsException {
	
		
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings)
			throws InvalidSettingsException {
		
		
	}

	@Override
	protected File getLockFile() {
		return new File(LOCKFILE);
	}
}


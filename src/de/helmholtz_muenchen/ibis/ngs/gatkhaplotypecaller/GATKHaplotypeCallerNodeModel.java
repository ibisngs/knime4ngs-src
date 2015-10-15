package de.helmholtz_muenchen.ibis.ngs.gatkhaplotypecaller;


import java.io.File;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.GVCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GATKHaplotypeCaller.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKHaplotypeCallerNodeModel extends GATKNodeModel {

	private String OUTFILE, LOCKFILE;
    private int bam_index;
	
	protected GATKHaplotypeCallerNodeModel(int INPORTS, int OUTPORTS) {
		super(OptionalPorts.createOPOs(INPORTS), OptionalPorts.createOPOs(OUTPORTS));
	}

	@Override
	protected String getCommandParameters(final BufferedDataTable[] inData) throws InvalidSettingsException {
		
		/**
    	 * Check INFILE
    	 */
    	String INFILE;
    	try{
    		INFILE = inData[0].iterator().next().getCell(bam_index).toString();
    		if(!INFILE.endsWith(".bam")){
    			throw new InvalidSettingsException("A cell of the input table has to be the path to BAM infile but it is "+INFILE);
    		}
    	}catch(IndexOutOfBoundsException e){
    			throw new InvalidSettingsException("A cell of the input table has to be the path to BAM infile but it is empty.");
    	}
		ArrayList<String> command 	= new ArrayList<String>();

		command.add("-I "+INFILE);
		command.add("--emitRefConfidence GVCF");
		command.add("--variant_index_type LINEAR");
		command.add("--variant_index_parameter 128000");
		
		this.OUTFILE = IO.replaceFileExtension(INFILE, ".gvcf");
		this.LOCKFILE = IO.replaceFileExtension(INFILE, SuccessfulRunChecker.LOCK_ENDING);
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

	@Override
	protected boolean checkInputCellType(DataTableSpec[] inSpecs) {
		bam_index = -1;
		
		for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
    		if(inSpecs[0].getColumnSpec(i).getType().toString().equals("BAMCell")) {
    			bam_index = i;
    		}
    	}
		return (bam_index>-1);
	}

	@Override
	protected DataType getOutColType() {
		return GVCFCell.TYPE;
	}
}


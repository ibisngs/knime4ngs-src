package de.helmholtz_muenchen.ibis.ngs.gatkcombinegvcfs;

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

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.GVCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of CombineGVCFs.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class CombineGVCFsNodeModel extends GATKNodeModel {
    
	private String OUTFILE, LOCKFILE;
	private int gvcf_index;
    
    /**
     * Constructor for the node model.
     */
    protected CombineGVCFsNodeModel() {
    
        super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
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
			String INFILE = row.getCell(gvcf_index).toString();
			
			if(first){
				this.OUTFILE = IO.replaceFileExtension(INFILE, ".ALLVARIANTS.gvcf");
				this.LOCKFILE = IO.replaceFileExtension(INFILE, SuccessfulRunChecker.LOCK_ENDING);
				first=false;
			}

			command.add("--variant "+INFILE);

		}

		return StringUtils.join(command, " ");
	}

	@Override
	protected String getCommandWalker() {
		return "CombineGVCFs";
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
		return new File(this.LOCKFILE);
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
		return GVCFCell.TYPE;
	}

}


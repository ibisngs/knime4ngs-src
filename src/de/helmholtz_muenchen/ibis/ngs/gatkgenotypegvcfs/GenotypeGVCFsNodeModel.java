package de.helmholtz_muenchen.ibis.ngs.gatkgenotypegvcfs;

import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataRow;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;

/**
 * This is the model implementation of GenotypeGVCFs.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GenotypeGVCFsNodeModel extends GATKNodeModel {
    
	private String OUTFILE; 

    protected GenotypeGVCFsNodeModel(int INPORTS, int OUTPORTS) {
		super(1, 1);
	}

	/**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
         // TODO: generated method stub
    }


	@Override
	protected String getCommandParameters(BufferedDataTable[] inData) {
		
		Iterator <DataRow> it = inData[0].iterator();
		ArrayList<String> command 	= new ArrayList<String>();
		boolean first = true;
		while(it.hasNext()){
			DataRow row = it.next();
			String INFILE = row.getCell(0).toString();
			
			if(first){
				this.OUTFILE = IO.replaceFileExtension(INFILE, ".GenotypedVariants.vcf");
				first=false;
			}

			command.add("--variant "+INFILE);

		}
		
		return StringUtils.join(command, " ");
	}

	@Override
	protected String getCommandWalker() {
		return "GenotypeGVCFs";
	}

	@Override
	protected String getOutfile() {
		return OUTFILE;
	}

	@Override
	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
			throws InvalidSettingsException {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings)
			throws InvalidSettingsException {
		// TODO Auto-generated method stub
		
	}
}


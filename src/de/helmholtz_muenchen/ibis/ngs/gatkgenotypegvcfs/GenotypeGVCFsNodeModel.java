package de.helmholtz_muenchen.ibis.ngs.gatkgenotypegvcfs;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataRow;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GenotypeGVCFs.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GenotypeGVCFsNodeModel extends GATKNodeModel {
    
	public static final String CFGKEY_NT_FILE = "NT";
//	private static final String REFERENCE = "Reference";

	private final SettingsModelIntegerBounded m_NT = new SettingsModelIntegerBounded(GenotypeGVCFsNodeModel.CFGKEY_NT_FILE, 1, 1, Integer.MAX_VALUE);
	
	private String OUTFILE; 

    protected GenotypeGVCFsNodeModel(int INPORTS, int OUTPORTS) {
		super(OptionalPorts.createOPOs(INPORTS), OptionalPorts.createOPOs(OUTPORTS));
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
		
//		String refGenome = getAvailableFlowVariables().get(REFERENCE).getStringValue();

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
//    	command.add("-R "+refGenome);
		command.add("-nt "+m_NT.getIntValue());
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
		m_NT.saveSettingsTo(settings);
		
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_NT.loadSettingsFrom(settings);
		
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_NT.validateSettings(settings);
		
	}
	@Override
	protected File getLockFile() {
		// TODO Auto-generated method stub
		return null;
	}
}


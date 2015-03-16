package de.helmholtz_muenchen.ibis.ngs.gatkcombinegvcfs;

import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataRow;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;

/**
 * This is the model implementation of CombineGVCFs.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class CombineGVCFsNodeModel extends GATKNodeModel {
    
	public static final String CFGKEY_BED_FILE 			= "BED_FILE";
	public static final String CFGKEY_BED_FILE_CHECKBOX = "BED_FILE_CHECKBOX";
	
	private final SettingsModelString m_BED_FILE = new SettingsModelString(CombineGVCFsNodeModel.CFGKEY_BED_FILE, "");
    private final SettingsModelBoolean m_BED_FILE_CHECKBOX = new SettingsModelBoolean(CombineGVCFsNodeModel.CFGKEY_BED_FILE_CHECKBOX, false);
	
	private String OUTFILE; 
    
    /**
     * Constructor for the node model.
     */
    protected CombineGVCFsNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(1, 1);
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
				this.OUTFILE = IO.replaceFileExtension(INFILE, ".ALLVARIANTS.vcf");
				first=false;
			}

			command.add("--variant "+INFILE);

		}
		if(m_BED_FILE_CHECKBOX.getBooleanValue()){
			command.add("-L "+m_BED_FILE.getStringValue());
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


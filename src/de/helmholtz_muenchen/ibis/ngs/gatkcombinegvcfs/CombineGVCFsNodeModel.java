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
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.CheckUtils;

import de.helmholtz_muenchen.ibis.utils.IO;
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
    
	static final String CFGKEY_OUTFOLDER = "outfolder";
    private final SettingsModelString m_OUTFOLDER = new SettingsModelString(CFGKEY_OUTFOLDER, "");
	
	
	private String OUTFILE;
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
				OUTFILE = m_OUTFOLDER.getStringValue()+ System.getProperty("file.separator")+ new File(INFILE).getName(); 
				OUTFILE = IO.replaceFileExtension(OUTFILE, ".ALLVARIANTS.gvcf");
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
	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
		m_OUTFOLDER.saveSettingsTo(settings);
	
	}

	@Override
	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
			throws InvalidSettingsException {
		m_OUTFOLDER.loadSettingsFrom(settings);
	
	}

	@Override
	protected void validateExtraSettings(NodeSettingsRO settings)
			throws InvalidSettingsException {
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
		return GVCFCell.TYPE;
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
	


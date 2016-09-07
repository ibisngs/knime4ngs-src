/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
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
        addSetting(m_OUTFOLDER);
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
	


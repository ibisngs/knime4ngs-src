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
package de.helmholtz_muenchen.ibis.ngs.gatkgenotypegvcfs;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GenotypeGVCFs.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GenotypeGVCFsNodeModel extends GATKNodeModel {
    
	static final String CFGKEY_NT_FILE = "NT";
	static final String CFGKEY_OUTFOLDER = "outfolder";

	private final SettingsModelIntegerBounded m_NT = new SettingsModelIntegerBounded(GenotypeGVCFsNodeModel.CFGKEY_NT_FILE, 1, 1, Integer.MAX_VALUE);
	private final SettingsModelString m_OUTFOLDER = new SettingsModelString(CFGKEY_OUTFOLDER, "");
	
	private String OUTFILE;
	private int gvcf_index;

    protected GenotypeGVCFsNodeModel(int INPORTS, int OUTPORTS) {
		super(OptionalPorts.createOPOs(INPORTS), OptionalPorts.createOPOs(OUTPORTS));
		addSetting(m_NT);
		addSetting(m_OUTFOLDER);
   	}

	@Override
	protected String getCommandParameters(BufferedDataTable[] inData) throws InvalidSettingsException {
		
		Iterator <DataRow> it = inData[0].iterator();
		ArrayList<String> command 	= new ArrayList<String>();
		boolean first = true;
		while(it.hasNext()){
			DataRow row = it.next();
			String INFILE = row.getCell(gvcf_index).toString();
			
			if(first){
				String outfolder = m_OUTFOLDER.getStringValue();
				if(outfolder.equals("") || outfolder == null) {
					outfolder = new File(INFILE).getParent();
				}
				if(!outfolder.endsWith(System.getProperty("file.separator"))) {
					outfolder = outfolder  + System.getProperty("file.separator");
				}
				OUTFILE = outfolder + new File(INFILE).getName();
				OUTFILE = IO.replaceFileExtension(OUTFILE, ".GenotypedVariants.vcf");
				first=false;
			}

			command.add("--variant "+INFILE);

		}

		command.add("-nt "+m_NT.getIntValue());
		return StringUtils.join(command, " ");
	}

	@Override
	protected String getCommandWalker() {
		return "GenotypeGVCFs";
	}


//	@Override
//	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
//		m_NT.saveSettingsTo(settings);
//		m_OUTFOLDER.saveSettingsTo(settings);
//		
//	}
//
//	@Override
//	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings)
//			throws InvalidSettingsException {
//		m_NT.loadSettingsFrom(settings);
//		m_OUTFOLDER.loadSettingsFrom(settings);
//		
//	}
//
//	@Override
//	protected void validateExtraSettings(NodeSettingsRO settings)
//			throws InvalidSettingsException {
//		m_NT.validateSettings(settings);
//		m_OUTFOLDER.validateSettings(settings);
//	}
	
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
		return VCFCell.TYPE;
	}

	@Override
	protected void extraConfig() throws InvalidSettingsException {
//		String outfolder_warning = CheckUtils.checkDestinationDirectory(m_OUTFOLDER.getStringValue());
//		if(outfolder_warning!=null) {
//			setWarningMessage(outfolder_warning);
//		}
	}
	
	@Override
	protected String getOutfile() {
		return OUTFILE;
	}
}


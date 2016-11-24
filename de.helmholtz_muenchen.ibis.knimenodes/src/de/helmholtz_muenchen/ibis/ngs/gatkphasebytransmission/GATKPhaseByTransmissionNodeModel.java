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
package de.helmholtz_muenchen.ibis.ngs.gatkphasebytransmission;


import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GATKPhaseByTransmission.
 * 
 *
 * @author Maximilian Hastreiter
 * @author Tim Jeske
 */
public class GATKPhaseByTransmissionNodeModel extends GATKNodeModel {

	/**
	 * Config Keys
	 */
	public static final String CFGKEY_PED_FILE = "PED";
	public static final String CFGKEY_DENOVOPRIOR = "deNovoPrior";

	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_PED_FILE = new SettingsModelString(CFGKEY_PED_FILE, ""); 
    private final SettingsModelString m_DENOVO_PRIOR = new SettingsModelString(CFGKEY_DENOVOPRIOR, "1.0E-8");
    
	//The Output Col Names
	public static final String OUT_COL1 = "PhasedVCF";
		
	private String OUTFILE;
	private int vcf_index;
	
    /**
     * Constructor for the node model.
     */
    protected GATKPhaseByTransmissionNodeModel() {
    	super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
    	addSetting(m_PED_FILE);
	   	addSetting(m_DENOVO_PRIOR);
    }

	@Override
	protected String getCommandParameters(BufferedDataTable[] inData)
			throws InvalidSettingsException {
		/**
    	 * Check INFILE
    	 */
    	String INFILE;
    	try{
    		INFILE = inData[0].iterator().next().getCell(vcf_index).toString();
    		if(!INFILE.endsWith(".vcf")){
    			throw new InvalidSettingsException("A cell of the input table has to be the path to VCF infile but it is "+INFILE);
    		}
    	}catch(IndexOutOfBoundsException e){
    			throw new InvalidSettingsException("A cell of the input table has to be the path to VCF infile but it is empty.");
    	}
    	
    	ArrayList<String> command = new ArrayList<String>();
    	command.add("-V "+INFILE);
    	
    	OUTFILE = IO.replaceFileExtension(INFILE, ".PhasedByTransmission.vcf");
    	
    	command.add("-prior "+m_DENOVO_PRIOR.getStringValue());
    	command.add("-ped "+m_PED_FILE.getStringValue());
    	
    	return StringUtils.join(command, " ");
	}
	

	@Override
	protected String getCommandWalker() {
		return "PhaseByTransmission";
	}


	@Override
	protected String getOutfile() {
		return OUTFILE;
	}

	@Override
	protected boolean checkInputCellType(DataTableSpec[] inSpecs) {
		vcf_index = -1;
		
		for(int i = 0; i < inSpecs[0].getNumColumns(); i++) {
    		if(inSpecs[0].getColumnSpec(i).getType().toString().equals("VCFCell")) {
    			vcf_index = i;
    		}
    	}
		return (vcf_index>-1);
	}

	@Override
	protected DataType getOutColType() {
		return VCFCell.TYPE;
	}

	@Override
	protected void extraConfig() throws InvalidSettingsException {
		// TODO Auto-generated method stub
		
	}
}


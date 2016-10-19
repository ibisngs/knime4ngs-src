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
package de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants;

import java.io.File;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GATKSelectVariants.
 * 
 *
 * @author Maximilian Hastreiter
 */
public abstract class SelectVariantsNodeModel extends GATKNodeModel {
    
	/**
	 * Config Keys
	 */
	public static final String CFGKEY_FILTERSTRING = "FILTERSTRING";
	
	/**
	 * The SettingsModels
	 */
	
    private final SettingsModelString m_FILTERSTRING = new SettingsModelString(CFGKEY_FILTERSTRING, "");
    
	private String OUTFILE, LOCKFILE;
	private int vcf_index;
	
	//The Output Col Names
	public static final String OUT_COL1_TABLE1 = "FilteredVCF";
    
    /**
     * Constructor for the node model.
     */
    protected SelectVariantsNodeModel(int INPORTS, int OUTPORTS) {
        super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
        addSetting(m_FILTERSTRING);
    }

    protected String getCommandParameters(final BufferedDataTable[] inData) throws InvalidSettingsException {
    	ArrayList<String> command = new ArrayList<String>();
    	
    	String INFILE = inData[0].iterator().next().getCell(vcf_index).toString();
    	try{
    		INFILE = inData[0].iterator().next().getCell(vcf_index).toString();
    		if(!INFILE.endsWith(".vcf")){
    			throw new InvalidSettingsException("A cell of the input table has to be the path to VCF infile but it is "+INFILE);
    		}
    	}catch(IndexOutOfBoundsException e){
    			throw new InvalidSettingsException("A cell of the input table has to be the path to VCF infile but it is empty.");
    	}
    	command.add("-V "+INFILE);
    	command.add(getCommandParameters());
    	this.OUTFILE = IO.replaceFileExtension(INFILE, getOutfileSuffix());
		this.LOCKFILE = IO.replaceFileExtension(INFILE, SuccessfulRunChecker.LOCK_ENDING);
    	
    	return StringUtils.join(command, " ");
    }
    
    protected File getLockFile() {
    	return new File(this.LOCKFILE);
    }
    
    
    protected String getOutfile() {
    	return this.OUTFILE;
    }

    protected String getCommandWalker() {
    	return "SelectVariants";
    }

//    protected void saveExtraSettingsTo(final NodeSettingsWO settings) {
//   	 m_FILTERSTRING.saveSettingsTo(settings);
//    }
//
//    protected void loadExtraValidatedSettingsFrom(final NodeSettingsRO settings)
//            throws InvalidSettingsException {
//       	 m_FILTERSTRING.loadSettingsFrom(settings);
//    }
//
//    protected void validateExtraSettings(final NodeSettingsRO settings)
//            throws InvalidSettingsException {
//       	 m_FILTERSTRING.validateSettings(settings);
//    }

    protected SettingsModelString getFILTERSTRINGModel(){
    	return m_FILTERSTRING;
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
    
    /****************************** ABSTRACT METHODS **********************************/
    /**
     * Provides the node specific filter settings
     * @return
     */
    protected abstract String getCommandParameters();
    
    /**
     * Provides the outfile suffix
     * @return
     */
    protected abstract String getOutfileSuffix();
    
}


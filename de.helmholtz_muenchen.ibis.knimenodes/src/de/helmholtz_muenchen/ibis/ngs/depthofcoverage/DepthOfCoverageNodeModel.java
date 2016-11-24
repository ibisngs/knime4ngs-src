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
package de.helmholtz_muenchen.ibis.ngs.depthofcoverage;


import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;



/**
 * This is the model implementation of DepthOfCoverage.
 * 
 *
 * @author Maximilian Hastreiter
 * @author Tim Jeske
 */
public class DepthOfCoverageNodeModel extends GATKNodeModel {
    
	 protected static final NodeLogger logger = NodeLogger.getLogger(DepthOfCoverageNodeModel.class);
	
	/**
	 * Config Keys
	 */
	public static final String CFGKEY_INFILE = "infile";
	public static final String CFGKEY_EXTRAFILTERS = "extrafilters";
	public static final String CFGKEY_FILESUFFIX = "suffix";
	
	/**
	 * Node Models
	 */
	private final SettingsModelOptionalString m_extrafilters = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_EXTRAFILTERS,"",false);
	private final SettingsModelString m_filesuffix = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_FILESUFFIX,"DoC");
	
	private String OUTFILE;
	private int bam_index;

    /**
     * Constructor for the node model.
     */
    protected DepthOfCoverageNodeModel() {
        super(OptionalPorts.createOPOs(1), OptionalPorts.createOPOs(1));
        addSetting(m_extrafilters);
        addSetting(m_filesuffix);
    }

    @Override
	protected String getCommandParameters(BufferedDataTable[] inData) throws InvalidSettingsException {

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
		
    	ArrayList<String> command = new ArrayList<String>();
		
		String fileSuffix = m_filesuffix.getStringValue();
		if(fileSuffix.equals("")) {
			fileSuffix = "DoC";
		}
		
		this.OUTFILE = IO.replaceFileExtension(INFILE, fileSuffix);
		command.add("-I "+INFILE);

		
		if(!m_extrafilters.getStringValue().equals("")){
        	String[] filters = m_extrafilters.getStringValue().split(",");
        	for(String filter : filters){
        		command.add(" --read_filter "+filter);
        	}
    	}
		
		DepthOfCoverageNodeModel.logger.info("Running GATK DepthOfCoverage...");
    	DepthOfCoverageNodeModel.logger.info("Log files can be found in "+OUTFILE+".out.log and "+OUTFILE+".err.log");
		
		return StringUtils.join(command, " ");
	}

	@Override
	protected String getCommandWalker() {
		return "DepthOfCoverage";
	}


	@Override
	protected String getOutfile() {
		return this.OUTFILE;
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
		return FileCell.TYPE;
	}

	@Override
	protected void extraConfig() throws InvalidSettingsException {
		// TODO Auto-generated method stub
		
	}
}


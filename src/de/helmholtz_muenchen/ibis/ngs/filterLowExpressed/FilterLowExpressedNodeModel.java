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
package de.helmholtz_muenchen.ibis.ngs.filterLowExpressed;

import java.io.File;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;

/**
 * This is the model implementation of FilterLowExpressed.
 * 
 *
 * @author Michael Kluge
 */
public class FilterLowExpressedNodeModel extends RNodeModel {

    // keys for SettingsModels
    protected static final String CFGKEY_KEEP_READS 	= "keepReads";
    protected static final String CFGKEY_KEEP_FRACTION 	= "keepFraction";
    protected static final String CFGKEY_BOTH_SEPERATE	= "bothConditionsSeperate";
    public static final String CFGKEY_MODE			= "operationModeFilter"; 
    
    // initial default values for SettingsModels  
    protected static final int DEFAULT_KEEP_READS = 10;
    protected static final double DEFAULT_KEEP_FRATION = 0.5;
    protected static final boolean DEFAULT_BOTH_SEP = true;
    public static final String DEFAULT_MODE = "Average over all samples";
    
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelInteger SET_KEEP_READS	= new SettingsModelInteger(CFGKEY_KEEP_READS, DEFAULT_KEEP_READS);
    private final SettingsModelDouble SET_KEEP_FRACTION	= new SettingsModelDouble(CFGKEY_KEEP_FRACTION, DEFAULT_KEEP_FRATION);
    private final SettingsModelBoolean SET_BOTH_SEP		= new SettingsModelBoolean(CFGKEY_BOTH_SEPERATE, DEFAULT_BOTH_SEP);
    private final SettingsModelString SET_MODE 			= new SettingsModelString(CFGKEY_MODE, DEFAULT_MODE);
    
    // the logger instance
	@SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(FilterLowExpressedNodeModel.class);
	
	private static final String SCRIPT_PATH = "ngs" + File.separatorChar + "de" + File.separatorChar + "filterLowExpressed.R";
	
    /**
     * Constructor for the node model.
     */
	protected FilterLowExpressedNodeModel() {
		super(2, 1, SCRIPT_PATH, new String[]{"--countTable", "--annotationFile"}, new String[]{"--output"});
		this.addSetting(SET_MODE);
		this.addSetting(SET_KEEP_READS);
		this.addSetting(SET_KEEP_FRACTION);
		this.addSetting(SET_BOTH_SEP);
		
    	// disable for standard mode
		SET_KEEP_FRACTION.setEnabled(false);
	}
	
    /**
     * {@inheritDoc}
     * @throws Exception 
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{    	
		BufferedDataTable[] out = super.execute(inData, exec);
		out[0] = exec.createSpecReplacerTable(out[0], this.getSpec(inData[0].getDataTableSpec())); // parse cell types

		return(out);
	}
    
    /**
     * get specs of table
     * @param inSpec
     * @return
     */
    private DataTableSpec getSpec(DataTableSpec inSpec) {
    	return inSpec;
	}
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	// add values set in GUI
    	this.addArgument("--keepReads", this.SET_KEEP_READS.getIntValue());
    	this.addArgument("--keepFraction", this.SET_KEEP_FRACTION.getDoubleValue());
    	this.addArgument("--bothConditionsSeperate", this.SET_BOTH_SEP.getBooleanValue() ? "1" : "0");
    	this.addArgument("--filterMode", this.SET_MODE.getStringValue().equals(DEFAULT_MODE) ? "1" : "0");
    	
    	return new DataTableSpec[]{getSpec(inSpecs[0])};
    }
}


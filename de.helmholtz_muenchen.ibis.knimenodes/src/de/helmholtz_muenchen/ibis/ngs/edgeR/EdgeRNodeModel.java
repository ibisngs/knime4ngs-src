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
package de.helmholtz_muenchen.ibis.ngs.edgeR;

import java.io.File;
import java.util.ArrayList;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.CompatibilityChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeModel;

/**
 * This is the model implementation of EdgeR.
 * 
 *
 * @author Michael Kluge
 */
public class EdgeRNodeModel extends RNodeModel {
	
    // keys for SettingsModels
    protected static final String CFGKEY_CORRECTION_METHOD 			= "PvalueCorrection";
    protected static final String CFGKEY_NORMALIZE_METHOD_FACTOR 	= "NormalizeMethodFactor";

    // initial default values for SettingsModels
    protected static final String DEFAULT_CORRECTION_METHOD = "BH";
    protected static final String DEFAULT_NORMALIZE_METHOD_FACTOR = "TMM";
    
	// definition of SettingsModel (all prefixed with SET)
    private final SettingsModelString SET_CORRECTION	= new SettingsModelString(CFGKEY_CORRECTION_METHOD, DEFAULT_CORRECTION_METHOD);
    private final SettingsModelString SET_NORM_FACTOR	= new SettingsModelString(CFGKEY_NORMALIZE_METHOD_FACTOR, DEFAULT_NORMALIZE_METHOD_FACTOR);
 
    // the logger instance
	@SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(EdgeRNodeModel.class);
	
	private static final String SCRIPT_PATH = "ngs" + File.separatorChar + "de" + File.separatorChar + "edgeR.R";
	
	// static vars for all availible methods
	protected static final ArrayList<String> CORRECTION_METHODS = new ArrayList<String>();
	protected static final ArrayList<String> NORM_FACTORS = new ArrayList<String>();
	
	static {
		CORRECTION_METHODS.add("BH");
		CORRECTION_METHODS.add("bonferroni");
		CORRECTION_METHODS.add("BY");
		CORRECTION_METHODS.add("fdr");
		CORRECTION_METHODS.add("hochberg");
		CORRECTION_METHODS.add("holm");
		CORRECTION_METHODS.add("hommel");
		CORRECTION_METHODS.add("none");
				
		NORM_FACTORS.add("RLE");
		NORM_FACTORS.add("TMM");
		NORM_FACTORS.add("upperquartile");
		NORM_FACTORS.add("none");
	}

    /**
     * Constructor for the node model.
     */
	protected EdgeRNodeModel() {
		super(2, 1, SCRIPT_PATH, new String[]{"--countTable", "--annotationFile"}, new String[]{"--output"}, 1);
		this.addSetting(SET_CORRECTION);
		this.addSetting(SET_NORM_FACTOR);
	}
	
    /**
     * {@inheritDoc}
     * @throws Exception 
     */
    @Override
	protected BufferedDataTable[] execute(final BufferedDataTable[] inData, final ExecutionContext exec) throws Exception{
    	
    	//Check input table integrity
    	CompatibilityChecker.inDataCheck(inData);
    	
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
    	DataColumnSpec[] specs = DataTableSpec.createColumnSpecs(new String[]{"ID", "log2FC", "aveLog2CPM", "PValue", "adj.PValue"}, new DataType[]{DataType.getType(StringCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class), DataType.getType(DoubleCell.class)});
    	return(new DataTableSpec(specs));
	}
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	
    	CompatibilityChecker.InputFileNoColumns(inSpecs);
    	
    	// add values set in GUI
    	this.addArgument("--normFactors", this.SET_NORM_FACTOR.getStringValue());
    	this.addArgument("--correctPvalue", this.SET_CORRECTION.getStringValue());

    	return new DataTableSpec[]{getSpec(inSpecs[0])};
    }
}


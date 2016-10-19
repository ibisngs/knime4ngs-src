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


package de.helmholtz_muenchen.ibis.ngs.gatkraw;


import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.ngs.OptionalPorts;

/**
 * This is the model implementation of GATKRaw.
 * GATK node without given command walker
 *
 * @author Tim Jeske
 */
public class GATKRawNodeModel extends GATKNodeModel {
    
	static final String CFGKEY_WALKER = "command_walker";
	static final String CFGKEY_OUT = "outfile";
	
	private final SettingsModelString m_walker = new SettingsModelString(GATKRawNodeModel.CFGKEY_WALKER,"");
	private final SettingsModelString m_out = new SettingsModelString(GATKRawNodeModel.CFGKEY_OUT,"");
	
    /**
     * Constructor for the node model.
     */
    protected GATKRawNodeModel(int in, int out) {
    
        super(OptionalPorts.createOPOs(in), OptionalPorts.createOPOs(out));
        addSetting(m_walker);
		addSetting(m_out);
    }

	@Override
	protected String getCommandParameters(BufferedDataTable[] inData) throws InvalidSettingsException {
		return "";
	}

	@Override
	protected String getCommandWalker() {
		return m_walker.getStringValue();
	}


	@Override
	protected String getOutfile() {
		return m_out.getStringValue();
	}

	@Override
	protected boolean checkInputCellType(DataTableSpec[] inSpecs) {
		return true;
	}

	@Override
	protected DataType getOutColType() {
		return FileCell.TYPE;
	}

	@Override
	protected void extraConfig() throws InvalidSettingsException {
		
	}

//	@Override
//	protected void saveExtraSettingsTo(NodeSettingsWO settings) {
//		m_walker.saveSettingsTo(settings);
//		m_out.saveSettingsTo(settings);
//	}
//
//	@Override
//	protected void loadExtraValidatedSettingsFrom(NodeSettingsRO settings) throws InvalidSettingsException {
//		m_walker.loadSettingsFrom(settings);
//		m_out.loadSettingsFrom(settings);
//	}
//
//	@Override
//	protected void validateExtraSettings(NodeSettingsRO settings) throws InvalidSettingsException {
//		m_walker.validateSettings(settings);
//		m_out.validateSettings(settings);
//	}
}


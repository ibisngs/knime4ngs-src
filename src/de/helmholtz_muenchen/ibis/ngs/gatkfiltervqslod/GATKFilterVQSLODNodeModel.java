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


package de.helmholtz_muenchen.ibis.ngs.gatkfiltervqslod;


import org.knime.core.node.InvalidSettingsException;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants.SelectVariantsNodeModel;

/**
 * This is the model implementation of GATKFilterVQSLOD.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKFilterVQSLODNodeModel extends SelectVariantsNodeModel {

	protected GATKFilterVQSLODNodeModel(int INPORTS, int OUTPORTS) {
		super(INPORTS, OUTPORTS);
	
	}

	@Override
	protected String getCommandParameters() {
		
		double Cutoff = Double.parseDouble(getFILTERSTRINGModel().getStringValue());
		
		return "--select_expressions VQSLOD>="+Cutoff;
	}

	@Override
	protected String getOutfileSuffix() {
		return "VQSLOD_Filtered.vcf";
	}

	@Override
	protected void extraConfig() throws InvalidSettingsException {
		// TODO Auto-generated method stub
		
	}
}


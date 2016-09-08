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


package de.helmholtz_muenchen.ibis.ngs.gatkexcludevariants;

import org.knime.core.node.InvalidSettingsException;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.SelectVariants.SelectVariantsNodeModel;

/**
 * This is the model implementation of GATKExcludeVariants.
 * 
 *
 * @author Tim Jeske
 */
public class GATKExcludeVariantsNodeModel extends SelectVariantsNodeModel {
    
    /**
     * Constructor for the node model.
     */
    protected GATKExcludeVariantsNodeModel() {
        super(1, 1);
        getFILTERSTRINGModel().setStringValue("--excludeFiltered --excludeNonVariants");
    }

	@Override
	protected String getCommandParameters() {
		return getFILTERSTRINGModel().getStringValue();
	}

	@Override
	protected String getOutfileSuffix() {
		return ".PASS_variants.vcf";
	}

	@Override
	protected void extraConfig() throws InvalidSettingsException {
		// TODO Auto-generated method stub
		
	}
}


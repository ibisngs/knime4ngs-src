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
package de.helmholtz_muenchen.ibis.ngs.fastqc;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;


/**
 * <code>NodeDialog</code> for the "FastQC" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class FastQCNodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the FastQC node.
     */
    protected FastQCNodeDialog() {
    	super(FastQCell.TYPE.getPreferredValueClass(), FastQCell.TYPE.getPreferredValueClass());
    }

	@Override
	public void addToolDialogComponents() {
		// TODO Auto-generated method stub
		
	}
}


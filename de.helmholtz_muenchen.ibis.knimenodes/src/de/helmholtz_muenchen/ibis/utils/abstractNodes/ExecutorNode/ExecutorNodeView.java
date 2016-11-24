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
package de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode;


import java.awt.Dimension;

import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode.ScriptNodeModel;


/**
 * <code>ExecutorNodeView</code> which shows STDOUT and STDERR in the view
 * 
 * @author Jonas Zierer
 * @author Michael Kluge
 *
 */
public abstract class ExecutorNodeView<T extends ExecutorNodeModel> extends NodeView<T> {

    // panel which actually shows the logs
	private final LogPanel PANEL_LOG; 
    
    /**
     * Creates a new view.
     * 
     * @param nodeModel the model class: {@link ScriptNodeModel}
     */
	protected ExecutorNodeView(final T nodeModel) {
        super(nodeModel);
        PANEL_LOG = new LogPanel(nodeModel.getSTDOUT(), nodeModel.getSTDERR());
        PANEL_LOG.setPreferredSize(new Dimension(800, 600));
        setComponent(PANEL_LOG);    
        setViewTitleSuffix("log");
        setShowNODATALabel(false); // always show view, never show no data
    }

    /** {@inheritDoc} */
    @Override
    protected void modelChanged() {
    	PANEL_LOG.updateView();
    }
    
	@Override
	protected void onClose() {
	
	}

	@Override
	protected void onOpen() {

	}
}

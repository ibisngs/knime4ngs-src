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

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "FilterLowExpressed" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FilterLowExpressedNodeFactory 
        extends NodeFactory<FilterLowExpressedNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FilterLowExpressedNodeModel createNodeModel() {
        return new FilterLowExpressedNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public HTENodeView<FilterLowExpressedNodeModel> createNodeView(final int viewIndex, final FilterLowExpressedNodeModel nodeModel) {
        return new HTENodeView<FilterLowExpressedNodeModel>(nodeModel);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new FilterLowExpressedNodeDialog();
    }

}


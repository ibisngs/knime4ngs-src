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


package de.helmholtz_muenchen.ibis.ngs.fastSam2Bam;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "FastSam2Bam" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FastSam2BamNodeView extends NodeView<FastSam2BamNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link FastSam2BamNodeModel})
     */
    protected FastSam2BamNodeView(final FastSam2BamNodeModel nodeModel) {
        super(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() { }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() { }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() { }
}


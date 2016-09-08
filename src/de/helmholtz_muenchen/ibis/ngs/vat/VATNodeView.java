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

package de.helmholtz_muenchen.ibis.ngs.vat;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "VAT" Node.
 * 
 *
 * @author Marie-Sophie Friedl
 */
public class VATNodeView extends NodeView<VATNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link VATNodeModel})
     */
    protected VATNodeView(final VATNodeModel nodeModel) {
        super(nodeModel);

        // instantiate the components of the view here.

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() {

        // retrieve the new model from your nodemodel and 
        // update the view.
        VATNodeModel nodeModel = 
            (VATNodeModel)getNodeModel();
        assert nodeModel != null;
        
        // be aware of a possibly not executed nodeModel! The data you retrieve
        // from your nodemodel could be null, emtpy, or invalid in any kind.
        
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() {
    
        // things to do when closing the view
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() {

        // things to do when opening the view
    }

}


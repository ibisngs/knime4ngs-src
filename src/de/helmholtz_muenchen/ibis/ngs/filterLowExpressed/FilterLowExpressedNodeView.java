package de.helmholtz_muenchen.ibis.ngs.filterLowExpressed;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "FilterLowExpressed" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FilterLowExpressedNodeView extends NodeView<FilterLowExpressedNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link FilterLowExpressedNodeModel})
     */
    protected FilterLowExpressedNodeView(final FilterLowExpressedNodeModel nodeModel) {
        super(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() {
    }

}


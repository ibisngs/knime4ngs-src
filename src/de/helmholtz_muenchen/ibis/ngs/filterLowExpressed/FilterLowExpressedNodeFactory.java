package de.helmholtz_muenchen.ibis.ngs.filterLowExpressed;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

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
    public NodeView<FilterLowExpressedNodeModel> createNodeView(final int viewIndex, final FilterLowExpressedNodeModel nodeModel) {
        return new FilterLowExpressedNodeView(nodeModel);
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


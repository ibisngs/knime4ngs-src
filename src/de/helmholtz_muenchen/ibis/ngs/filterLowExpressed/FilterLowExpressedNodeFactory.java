package de.helmholtz_muenchen.ibis.ngs.filterLowExpressed;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

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
    public RNodeView<FilterLowExpressedNodeModel> createNodeView(final int viewIndex, final FilterLowExpressedNodeModel nodeModel) {
        return new RNodeView<FilterLowExpressedNodeModel>(nodeModel);
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


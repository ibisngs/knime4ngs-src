package de.helmholtz_muenchen.ibis.ngs.filtersummary;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "FilterSummary" Node.
 * 
 *
 * @author 
 */
public class FilterSummaryNodeFactory 
        extends NodeFactory<FilterSummaryNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FilterSummaryNodeModel createNodeModel() {
        return new FilterSummaryNodeModel();
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
    public NodeView<FilterSummaryNodeModel> createNodeView(final int viewIndex,
            final FilterSummaryNodeModel nodeModel) {
        return new FilterSummaryNodeView(nodeModel);
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
        return new FilterSummaryNodeDialog();
    }

}


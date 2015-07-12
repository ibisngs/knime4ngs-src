package de.helmholtz_muenchen.ibis.ngs.extendgenesummary;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "ExtendGeneSummary" Node.
 * 
 *
 * @author tim.jeske
 */
public class ExtendGeneSummaryNodeFactory 
        extends NodeFactory<ExtendGeneSummaryNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ExtendGeneSummaryNodeModel createNodeModel() {
        return new ExtendGeneSummaryNodeModel();
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
    public NodeView<ExtendGeneSummaryNodeModel> createNodeView(final int viewIndex,
            final ExtendGeneSummaryNodeModel nodeModel) {
        return new ExtendGeneSummaryNodeView(nodeModel);
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
        return new ExtendGeneSummaryNodeDialog();
    }

}


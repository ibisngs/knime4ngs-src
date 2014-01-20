package de.helmholtz_muenchen.ibis.ngs.submittoqueue;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "SubmitToQueue" Node.
 * 
 *
 * @author 
 */
public class SubmitToQueueNodeFactory 
        extends NodeFactory<SubmitToQueueNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SubmitToQueueNodeModel createNodeModel() {
        return new SubmitToQueueNodeModel();
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
    public NodeView<SubmitToQueueNodeModel> createNodeView(final int viewIndex,
            final SubmitToQueueNodeModel nodeModel) {
        return new SubmitToQueueNodeView(nodeModel);
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
        return new SubmitToQueueNodeDialog();
    }

}


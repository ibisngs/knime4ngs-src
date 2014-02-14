package de.helmholtz_muenchen.ibis.sleepnode;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "SleepNode" Node.
 * Node sleeps for the given time
 *
 * @author Jonas Zierer
 */
public class SleepNodeNodeFactory 
        extends NodeFactory<SleepNodeNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SleepNodeNodeModel createNodeModel() {
        return new SleepNodeNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<SleepNodeNodeModel> createNodeView(final int viewIndex, final SleepNodeNodeModel nodeModel) {
        return null;
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
        return new SleepNodeNodeDialog();
    }

}


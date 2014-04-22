package de.helmholtz_muenchen.ibis.misc.shuffleData;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "ShuffleData" Node.
 * Shuffle data within columns of the input data matrix.
 *
 * @author Jonas Zierer
 */
public class ShuffleDataNodeFactory extends NodeFactory<ShuffleDataNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ShuffleDataNodeModel createNodeModel() {
        return new ShuffleDataNodeModel();
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
    public NodeView<ShuffleDataNodeModel> createNodeView(final int viewIndex, final ShuffleDataNodeModel nodeModel) {
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
        return new ShuffleDataNodeDialog();
    }

}


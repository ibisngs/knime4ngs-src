package de.helmholtz_muenchen.ibis.knimenodesstarter;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "KnimeNodesStarter" Node.
 * 
 *
 * @author hastreiter
 */
public class KnimeNodesStarterNodeFactory 
        extends NodeFactory<KnimeNodesStarterNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public KnimeNodesStarterNodeModel createNodeModel() {
        return new KnimeNodesStarterNodeModel();
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
    public NodeView<KnimeNodesStarterNodeModel> createNodeView(final int viewIndex,
            final KnimeNodesStarterNodeModel nodeModel) {
        return new KnimeNodesStarterNodeView(nodeModel);
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
        return new KnimeNodesStarterNodeDialog();
    }

}


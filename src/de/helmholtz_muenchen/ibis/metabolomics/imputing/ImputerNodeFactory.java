package de.helmholtz_muenchen.ibis.metabolomics.imputing;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Imputer" Node.
 * 
 *
 * @author Jonas Zierer
 */
public class ImputerNodeFactory extends NodeFactory<ImputerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ImputerNodeModel createNodeModel() {
        return new ImputerNodeModel();
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
    public NodeView<ImputerNodeModel> createNodeView(final int viewIndex,
            final ImputerNodeModel nodeModel) {
        return new ImputerNodeView(nodeModel);
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
        return new ImputerNodeDialog();
    }

}


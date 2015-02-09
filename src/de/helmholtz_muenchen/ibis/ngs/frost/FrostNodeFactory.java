package de.helmholtz_muenchen.ibis.ngs.frost;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Frost" Node.
 * 
 *
 * @author Syeda Tanzeem Haque
 */
public class FrostNodeFactory 
        extends NodeFactory<FrostNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FrostNodeModel createNodeModel() {
        return new FrostNodeModel();
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
    public NodeView<FrostNodeModel> createNodeView(final int viewIndex,
            final FrostNodeModel nodeModel) {
        return new FrostNodeView(nodeModel);
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
        return new FrostNodeDialog();
    }

}


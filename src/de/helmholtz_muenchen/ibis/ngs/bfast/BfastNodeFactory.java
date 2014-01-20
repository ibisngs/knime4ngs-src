package de.helmholtz_muenchen.ibis.ngs.bfast;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Bfast" Node.
 * 
 *
 * @author 
 */
public class BfastNodeFactory 
        extends NodeFactory<BfastNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public BfastNodeModel createNodeModel() {
        return new BfastNodeModel();
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
    public NodeView<BfastNodeModel> createNodeView(final int viewIndex,
            final BfastNodeModel nodeModel) {
        return new BfastNodeView(nodeModel);
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
        return new BfastNodeDialog();
    }

}


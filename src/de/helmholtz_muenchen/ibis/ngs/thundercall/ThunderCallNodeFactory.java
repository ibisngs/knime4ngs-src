package de.helmholtz_muenchen.ibis.ngs.thundercall;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "ThunderCall" Node.
 * 
 *
 * @author Tanzeem Haque
 */
public class ThunderCallNodeFactory 
        extends NodeFactory<ThunderCallNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ThunderCallNodeModel createNodeModel() {
        return new ThunderCallNodeModel();
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
    public NodeView<ThunderCallNodeModel> createNodeView(final int viewIndex,
            final ThunderCallNodeModel nodeModel) {
        return new ThunderCallNodeView(nodeModel);
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
        return new ThunderCallNodeDialog();
    }

}


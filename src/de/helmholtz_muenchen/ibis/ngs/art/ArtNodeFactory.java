package de.helmholtz_muenchen.ibis.ngs.art;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Art" Node.
 * 
 *
 * @author Syeda Tanzeem Haque
 */
public class ArtNodeFactory 
        extends NodeFactory<ArtNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ArtNodeModel createNodeModel() {
        return new ArtNodeModel();
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
    public NodeView<ArtNodeModel> createNodeView(final int viewIndex,
            final ArtNodeModel nodeModel) {
        return new ArtNodeView(nodeModel);
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
        return new ArtNodeDialog();
    }

}


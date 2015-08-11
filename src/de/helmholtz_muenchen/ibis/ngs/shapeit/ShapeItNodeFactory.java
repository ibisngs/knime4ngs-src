package de.helmholtz_muenchen.ibis.ngs.shapeit;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "ShapeIt" Node.
 * 
 *
 * @author Tanzeem Haque
 */
public class ShapeItNodeFactory 
        extends NodeFactory<ShapeItNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ShapeItNodeModel createNodeModel() {
        return new ShapeItNodeModel();
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
    public NodeView<ShapeItNodeModel> createNodeView(final int viewIndex,
            final ShapeItNodeModel nodeModel) {
        return new ShapeItNodeView(nodeModel);
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
        return new ShapeItNodeDialog();
    }

}


package de.helmholtz_muenchen.ibis.ngs.edgeR;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "EdgeR" Node.
 * 
 *
 * @author Michael Kluge
 */
public class EdgeRNodeFactory 
        extends NodeFactory<EdgeRNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public EdgeRNodeModel createNodeModel() {
        return new EdgeRNodeModel();
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
    public NodeView<EdgeRNodeModel> createNodeView(final int viewIndex,
            final EdgeRNodeModel nodeModel) {
        return new EdgeRNodeView(nodeModel);
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
        return new EdgeRNodeDialog();
    }

}


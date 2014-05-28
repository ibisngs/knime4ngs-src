package de.helmholtz_muenchen.ibis.plotting.boxplot;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeFactory</code> for the "Boxplot" Node.
 * 
 *
 * @author Jonas Zierer
 */
public class BoxplotNodeFactory 
        extends NodeFactory<BoxplotNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public BoxplotNodeModel createNodeModel() {
        return new BoxplotNodeModel();
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
    public NodeView<BoxplotNodeModel> createNodeView(final int viewIndex, final BoxplotNodeModel nodeModel) {
        return new RNodeView<BoxplotNodeModel>(nodeModel);
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
        return new BoxplotNodeDialog();
    }

}


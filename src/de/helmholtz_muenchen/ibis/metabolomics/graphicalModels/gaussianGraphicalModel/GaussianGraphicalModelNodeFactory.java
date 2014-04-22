package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.gaussianGraphicalModel;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeFactory</code> for the "RNodeTemplate" Node.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class GaussianGraphicalModelNodeFactory extends NodeFactory<GaussianGraphicalModelNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GaussianGraphicalModelNodeModel createNodeModel() {
        return new GaussianGraphicalModelNodeModel();
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
    public NodeView<GaussianGraphicalModelNodeModel> createNodeView(final int viewIndex, final GaussianGraphicalModelNodeModel nodeModel) {
        return new RNodeView<GaussianGraphicalModelNodeModel>(nodeModel);
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
        return new GaussianGraphicalModelNodeDialog();
    }

}


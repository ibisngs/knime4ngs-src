package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.edgeRanking;

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
public class GraphicalModelsNodeFactory extends NodeFactory<GraphicalModelsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GraphicalModelsNodeModel createNodeModel() {
        return new GraphicalModelsNodeModel();
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
    public NodeView<GraphicalModelsNodeModel> createNodeView(final int viewIndex, final GraphicalModelsNodeModel nodeModel) {
        return new RNodeView<GraphicalModelsNodeModel>(nodeModel);
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
        return new GraphicalModelsNodeDialog();
    }

}


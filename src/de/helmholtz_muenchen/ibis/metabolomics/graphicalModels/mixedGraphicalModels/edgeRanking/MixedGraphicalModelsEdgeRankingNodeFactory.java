package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.edgeRanking;

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
public class MixedGraphicalModelsEdgeRankingNodeFactory extends NodeFactory<MixedGraphicalModelsEdgeRankingNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public MixedGraphicalModelsEdgeRankingNodeModel createNodeModel() {
        return new MixedGraphicalModelsEdgeRankingNodeModel();
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
    public NodeView<MixedGraphicalModelsEdgeRankingNodeModel> createNodeView(final int viewIndex, final MixedGraphicalModelsEdgeRankingNodeModel nodeModel) {
        return new RNodeView<MixedGraphicalModelsEdgeRankingNodeModel>(nodeModel);
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
        return new MixedGraphicalModelsEdgeRankingNodeDialog();
    }

}


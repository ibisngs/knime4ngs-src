package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.edgeRankingFDRbackground;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "RNodeTemplate" Node.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class MixedGraphicalModelsEdgeRankingBackgroundNodeFactory extends NodeFactory<MixedGraphicalModelsEdgeRankingBackgroundNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public MixedGraphicalModelsEdgeRankingBackgroundNodeModel createNodeModel() {
        return new MixedGraphicalModelsEdgeRankingBackgroundNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<MixedGraphicalModelsEdgeRankingBackgroundNodeModel> createNodeView(final int viewIndex, final MixedGraphicalModelsEdgeRankingBackgroundNodeModel nodeModel) {
        return null;
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
        return new MixedGraphicalModelsEdgeRankingBackgroundNodeDialog();
    }

}


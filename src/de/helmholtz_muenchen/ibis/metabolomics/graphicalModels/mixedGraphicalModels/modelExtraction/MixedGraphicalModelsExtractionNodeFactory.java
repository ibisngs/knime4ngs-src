package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.modelExtraction;

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
public class MixedGraphicalModelsExtractionNodeFactory extends NodeFactory<MixedGraphicalModelsExtractionNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public MixedGraphicalModelsExtractionNodeModel createNodeModel() {
        return new MixedGraphicalModelsExtractionNodeModel();
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
    public NodeView<MixedGraphicalModelsExtractionNodeModel> createNodeView(final int viewIndex, final MixedGraphicalModelsExtractionNodeModel nodeModel) {
        return new RNodeView<MixedGraphicalModelsExtractionNodeModel>(nodeModel);
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
        return new MixedGraphicalModelsExtractionNodeDialog();
    }

}


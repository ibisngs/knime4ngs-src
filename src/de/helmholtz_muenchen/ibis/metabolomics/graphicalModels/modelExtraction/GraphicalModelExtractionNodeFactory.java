package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.modelExtraction;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "RNodeTemplate" Node.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class GraphicalModelExtractionNodeFactory 
        extends NodeFactory<GraphicalModelExtractionNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GraphicalModelExtractionNodeModel createNodeModel() {
        return new GraphicalModelExtractionNodeModel();
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
    public NodeView<GraphicalModelExtractionNodeModel> createNodeView(final int viewIndex,
            final GraphicalModelExtractionNodeModel nodeModel) {
        return new GraphicalModelExtractionNodeView(nodeModel);
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
        return new GraphicalModelExtractionNodeDialog();
    }

}


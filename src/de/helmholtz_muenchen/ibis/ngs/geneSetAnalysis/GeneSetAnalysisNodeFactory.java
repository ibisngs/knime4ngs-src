package de.helmholtz_muenchen.ibis.ngs.geneSetAnalysis;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GeneSetAnalysis" Node.
 * 
 *
 * @author Tim Jeske
 */
public class GeneSetAnalysisNodeFactory 
        extends NodeFactory<GeneSetAnalysisNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GeneSetAnalysisNodeModel createNodeModel() {
        return new GeneSetAnalysisNodeModel();
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
    public NodeView<GeneSetAnalysisNodeModel> createNodeView(final int viewIndex,
            final GeneSetAnalysisNodeModel nodeModel) {
        return new GeneSetAnalysisNodeView(nodeModel);
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
        return new GeneSetAnalysisNodeDialog();
    }

}


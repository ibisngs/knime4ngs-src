package de.helmholtz_muenchen.ibis.ngs.geneBasedAnalysis;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GeneBasedAnalysis" Node.
 * 
 *
 * @author Tim Jeske
 */
public class GeneBasedAnalysisNodeFactory 
        extends NodeFactory<GeneBasedAnalysisNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GeneBasedAnalysisNodeModel createNodeModel() {
        return new GeneBasedAnalysisNodeModel();
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
    public NodeView<GeneBasedAnalysisNodeModel> createNodeView(final int viewIndex,
            final GeneBasedAnalysisNodeModel nodeModel) {
        return new GeneBasedAnalysisNodeView(nodeModel);
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
        return new GeneBasedAnalysisNodeDialog();
    }

}


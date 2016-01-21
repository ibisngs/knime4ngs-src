package de.helmholtz_muenchen.ibis.ngs.transcriptBasedAnalysis;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "TranscriptBasedAnalysis" Node.
 * 
 *
 * @author Tim Jeske
 */
public class TranscriptBasedAnalysisNodeFactory 
        extends NodeFactory<TranscriptBasedAnalysisNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public TranscriptBasedAnalysisNodeModel createNodeModel() {
        return new TranscriptBasedAnalysisNodeModel();
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
    public NodeView<TranscriptBasedAnalysisNodeModel> createNodeView(final int viewIndex,
            final TranscriptBasedAnalysisNodeModel nodeModel) {
        return new TranscriptBasedAnalysisNodeView(nodeModel);
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
        return new TranscriptBasedAnalysisNodeDialog();
    }

}


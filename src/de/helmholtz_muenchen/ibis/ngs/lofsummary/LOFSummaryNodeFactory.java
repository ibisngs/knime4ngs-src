package de.helmholtz_muenchen.ibis.ngs.lofsummary;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "LOFStatistics" Node.
 * 
 *
 * @author tim.jeske
 */
public class LOFSummaryNodeFactory 
        extends NodeFactory<LOFSummaryNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public LOFSummaryNodeModel createNodeModel() {
        return new LOFSummaryNodeModel();
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
    public NodeView<LOFSummaryNodeModel> createNodeView(final int viewIndex,
            final LOFSummaryNodeModel nodeModel) {
        return new LOFSummaryNodeView(nodeModel);
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
        return new LOFSummaryNodeDialog();
    }

}


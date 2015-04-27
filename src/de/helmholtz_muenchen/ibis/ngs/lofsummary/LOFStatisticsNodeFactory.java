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
public class LOFStatisticsNodeFactory 
        extends NodeFactory<LOFStatisticsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public LOFStatisticsNodeModel createNodeModel() {
        return new LOFStatisticsNodeModel();
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
    public NodeView<LOFStatisticsNodeModel> createNodeView(final int viewIndex,
            final LOFStatisticsNodeModel nodeModel) {
        return new LOFStatisticsNodeView(nodeModel);
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


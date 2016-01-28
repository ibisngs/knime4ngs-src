package de.helmholtz_muenchen.ibis.ngs.vepsummary;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "LOFStatistics" Node.
 * 
 *
 * @author tim.jeske
 */
public class VEPSummaryNodeFactory 
        extends NodeFactory<VEPSummaryNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public VEPSummaryNodeModel createNodeModel() {
        return new VEPSummaryNodeModel();
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
    public NodeView<VEPSummaryNodeModel> createNodeView(final int viewIndex,
            final VEPSummaryNodeModel nodeModel) {
        return new VEPSummaryNodeView(nodeModel);
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
        return new VEPSummaryNodeDialog();
    }

}


package de.helmholtz_muenchen.ibis.ngs.gatkphasebytransmission;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GATKPhaseByTransmission" Node.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKPhaseByTransmissionNodeFactory 
        extends NodeFactory<GATKPhaseByTransmissionNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKPhaseByTransmissionNodeModel createNodeModel() {
        return new GATKPhaseByTransmissionNodeModel();
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
    public NodeView<GATKPhaseByTransmissionNodeModel> createNodeView(final int viewIndex,
            final GATKPhaseByTransmissionNodeModel nodeModel) {
        return new GATKPhaseByTransmissionNodeView(nodeModel);
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
        return new GATKPhaseByTransmissionNodeDialog();
    }

}


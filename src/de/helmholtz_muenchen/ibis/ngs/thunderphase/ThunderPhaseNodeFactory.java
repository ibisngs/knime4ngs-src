package de.helmholtz_muenchen.ibis.ngs.thunderphase;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "ThunderPhase" Node.
 * 
 *
 * @author Tanzeem Haque
 */
public class ThunderPhaseNodeFactory 
        extends NodeFactory<ThunderPhaseNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ThunderPhaseNodeModel createNodeModel() {
        return new ThunderPhaseNodeModel();
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
    public NodeView<ThunderPhaseNodeModel> createNodeView(final int viewIndex,
            final ThunderPhaseNodeModel nodeModel) {
        return new ThunderPhaseNodeView(nodeModel);
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
        return new ThunderPhaseNodeDialog();
    }

}


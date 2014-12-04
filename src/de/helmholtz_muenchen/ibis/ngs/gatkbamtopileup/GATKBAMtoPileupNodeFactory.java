package de.helmholtz_muenchen.ibis.ngs.gatkbamtopileup;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GATKBAMtoPileup" Node.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKBAMtoPileupNodeFactory 
        extends NodeFactory<GATKBAMtoPileupNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKBAMtoPileupNodeModel createNodeModel() {
        return new GATKBAMtoPileupNodeModel();
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
    public NodeView<GATKBAMtoPileupNodeModel> createNodeView(final int viewIndex,
            final GATKBAMtoPileupNodeModel nodeModel) {
        return new GATKBAMtoPileupNodeView(nodeModel);
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
        return new GATKBAMtoPileupNodeDialog();
    }

}


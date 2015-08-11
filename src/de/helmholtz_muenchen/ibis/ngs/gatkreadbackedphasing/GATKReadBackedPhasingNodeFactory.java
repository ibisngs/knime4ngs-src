package de.helmholtz_muenchen.ibis.ngs.gatkreadbackedphasing;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GATKReadBackedPhasing" Node.
 * 
 *
 * @author Tanzeem Haque
 */
public class GATKReadBackedPhasingNodeFactory 
        extends NodeFactory<GATKReadBackedPhasingNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKReadBackedPhasingNodeModel createNodeModel() {
        return new GATKReadBackedPhasingNodeModel();
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
    public NodeView<GATKReadBackedPhasingNodeModel> createNodeView(final int viewIndex,
            final GATKReadBackedPhasingNodeModel nodeModel) {
        return new GATKReadBackedPhasingNodeView(nodeModel);
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
        return new GATKReadBackedPhasingNodeDialog();
    }

}


package de.helmholtz_muenchen.ibis.ngs.gatkunifiedgenotyper;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GATKUnifiedGenotyper" Node.
 * 
 *
 * @author 
 */
public class GATKUnifiedGenotyperNodeFactory 
        extends NodeFactory<GATKUnifiedGenotyperNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKUnifiedGenotyperNodeModel createNodeModel() {
        return new GATKUnifiedGenotyperNodeModel();
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
    public NodeView<GATKUnifiedGenotyperNodeModel> createNodeView(final int viewIndex,
            final GATKUnifiedGenotyperNodeModel nodeModel) {
        return new GATKUnifiedGenotyperNodeView(nodeModel);
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
        return new GATKUnifiedGenotyperNodeDialog();
    }

}


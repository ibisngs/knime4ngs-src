package de.helmholtz_muenchen.ibis.ngs.reordervcf;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "ReorderVCF" Node.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class ReorderVCFNodeFactory 
        extends NodeFactory<ReorderVCFNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ReorderVCFNodeModel createNodeModel() {
        return new ReorderVCFNodeModel();
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
    public NodeView<ReorderVCFNodeModel> createNodeView(final int viewIndex,
            final ReorderVCFNodeModel nodeModel) {
        return new ReorderVCFNodeView(nodeModel);
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
        return new ReorderVCFNodeDialog();
    }

}


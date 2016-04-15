package de.helmholtz_muenchen.ibis.ngs.gatkrealignment;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "GATKRealignment" Node.
 * 
 *
 * @author 
 */
public class GATKRealignmentNodeFactory 
        extends NodeFactory<GATKRealignmentNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKRealignmentNodeModel createNodeModel() {
        return new GATKRealignmentNodeModel();
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
    public NodeView<GATKRealignmentNodeModel> createNodeView(final int viewIndex,
            final GATKRealignmentNodeModel nodeModel) {
        return new HTENodeView<GATKRealignmentNodeModel>(nodeModel);
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
        return new GATKRealignmentNodeDialog();
    }

}


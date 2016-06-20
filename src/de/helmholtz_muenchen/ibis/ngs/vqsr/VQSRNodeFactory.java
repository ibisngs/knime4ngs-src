package de.helmholtz_muenchen.ibis.ngs.vqsr;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "VQSR" Node.
 * 
 *
 * @author 
 */
public class VQSRNodeFactory 
        extends NodeFactory<VQSRNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public VQSRNodeModel createNodeModel() {
        return new VQSRNodeModel();
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
    public NodeView<VQSRNodeModel> createNodeView(final int viewIndex,
            final VQSRNodeModel nodeModel) {
        return new HTENodeView<VQSRNodeModel>(nodeModel);
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
        return new VQSRNodeDialog();
    }

}


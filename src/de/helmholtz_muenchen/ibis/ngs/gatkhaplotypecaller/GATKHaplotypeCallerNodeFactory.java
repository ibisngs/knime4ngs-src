package de.helmholtz_muenchen.ibis.ngs.gatkhaplotypecaller;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "GATKHaplotypeCaller" Node.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKHaplotypeCallerNodeFactory 
        extends NodeFactory<GATKHaplotypeCallerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKHaplotypeCallerNodeModel createNodeModel() {
        return new GATKHaplotypeCallerNodeModel(1,1);
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
    public NodeView<GATKHaplotypeCallerNodeModel> createNodeView(final int viewIndex,
            final GATKHaplotypeCallerNodeModel nodeModel) {
        return new HTENodeView<GATKHaplotypeCallerNodeModel>(nodeModel);
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
        return new GATKHaplotypeCallerNodeDialog();
    }

}


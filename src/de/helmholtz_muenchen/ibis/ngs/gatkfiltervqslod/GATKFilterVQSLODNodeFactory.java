package de.helmholtz_muenchen.ibis.ngs.gatkfiltervqslod;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GATKFilterVQSLOD" Node.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKFilterVQSLODNodeFactory 
        extends NodeFactory<GATKFilterVQSLODNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKFilterVQSLODNodeModel createNodeModel() {
        return new GATKFilterVQSLODNodeModel(1,1);
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
    public NodeView<GATKFilterVQSLODNodeModel> createNodeView(final int viewIndex,
            final GATKFilterVQSLODNodeModel nodeModel) {
        return new GATKFilterVQSLODNodeView(nodeModel);
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
        return new GATKFilterVQSLODNodeDialog();
    }

}


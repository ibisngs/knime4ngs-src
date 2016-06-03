package de.helmholtz_muenchen.ibis.ngs.gatkvariantfiltration;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "GATKVariantFiltration" Node.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GATKVariantFiltrationNodeFactory 
        extends NodeFactory<GATKVariantFiltrationNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKVariantFiltrationNodeModel createNodeModel() {
        return new GATKVariantFiltrationNodeModel();
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
    public NodeView<GATKVariantFiltrationNodeModel> createNodeView(final int viewIndex,
            final GATKVariantFiltrationNodeModel nodeModel) {
        return new HTENodeView<GATKVariantFiltrationNodeModel>(nodeModel);
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
        return new GATKVariantFiltrationNodeDialog();
    }

}


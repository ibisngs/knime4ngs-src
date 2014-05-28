package de.helmholtz_muenchen.ibis.misc.outlierRemoval;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeFactory</code> for the "OutlierRemoval" Node.
 * Remove Outliers from Data Columns
 *
 * @author Jonas Zierer
 */
public class OutlierRemovalNodeFactory  extends NodeFactory<OutlierRemovalNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public OutlierRemovalNodeModel createNodeModel() {
        return new OutlierRemovalNodeModel();
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
    public NodeView<OutlierRemovalNodeModel> createNodeView(final int viewIndex, final OutlierRemovalNodeModel nodeModel) {
        return new RNodeView<OutlierRemovalNodeModel>(nodeModel);
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
        return new OutlierRemovalNodeDialog();
    }

}


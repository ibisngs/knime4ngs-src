package de.helmholtz_muenchen.ibis.misc.pca;

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
public class PCANodeFactory  extends NodeFactory<PCANodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public PCANodeModel createNodeModel() {
        return new PCANodeModel();
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
    public NodeView<PCANodeModel> createNodeView(final int viewIndex, final PCANodeModel nodeModel) {
        return new RNodeView<PCANodeModel>(nodeModel);
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
        return new PCANodeDialog();
    }

}


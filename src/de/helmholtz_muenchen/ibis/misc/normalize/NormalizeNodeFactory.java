package de.helmholtz_muenchen.ibis.misc.normalize;

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
public class NormalizeNodeFactory  extends NodeFactory<NormalizeNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public NormalizeNodeModel createNodeModel() {
        return new NormalizeNodeModel();
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
    public NodeView<NormalizeNodeModel> createNodeView(final int viewIndex, final NormalizeNodeModel nodeModel) {
        return new RNodeView<NormalizeNodeModel>(nodeModel);
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
        return new NormalizeNodeDialog();
    }

}


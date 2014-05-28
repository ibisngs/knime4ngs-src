package de.helmholtz_muenchen.ibis.misc.regression.batch;

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
public class RegressionBatchNodeFactory  extends NodeFactory<RegressionBatchNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public RegressionBatchNodeModel createNodeModel() {
        return new RegressionBatchNodeModel();
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
    public NodeView<RegressionBatchNodeModel> createNodeView(final int viewIndex, final RegressionBatchNodeModel nodeModel) {
        return new RNodeView<RegressionBatchNodeModel>(nodeModel);
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
        return new RegressionBatchNodeDialog();
    }

}

